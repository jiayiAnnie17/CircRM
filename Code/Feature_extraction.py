from __future__ import absolute_import
import os
import sys
import re
import h5py
from statsmodels import robust
import numpy as np
import csv
import pandas as pd
import argparse
import multiprocessing
from tqdm import tqdm
from itertools import chain

def extract_read_signals(file_dir):
    """Extract alignment and corrected signal from a single fast5 file"""
    read_id = file_dir.split('/')[-1].replace('.fast5', '')
    fast5_data = h5py.File(file_dir, 'r')
    group = FLAGS.group

    try:
        corr_data = fast5_data['/Analyses/'+group+'/BaseCalled_template/Events'][()]
    except:
        raise RuntimeError('Corrected data not found.')

    aln_data = fast5_data['/Analyses/'+ group + '/BaseCalled_template/Alignment']
    aln_data = dict(list(aln_data.attrs.items()))
    fast5_data.close()
    return (read_id, aln_data, corr_data)


def onehot_encode(sequence):
    """One-hot encode a RNA sequence"""
    encoding = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'C': [0, 0, 1, 0], 'G': [0, 0, 0, 1]}
    try:
        encoded= [encoding[base] for base in sequence]
        flat_encoded = list(chain.from_iterable(encoded))
        return [flat_encoded[i:i + 20] for i in range(0, len(flat_encoded), 20)]
    except KeyError as e:
        print(f"Unexpected base in sequence: {e}")
        return None


def extract_kmer_signals(read_id, aln_data, corr_data, kmer_filter):
    """Extract features for all 5-mers matching the filter"""
    event_mean = corr_data['norm_mean']
    event_stdev = corr_data['norm_stdev']
    event_lengths = corr_data['length']
    event_bases = corr_data['base']

    seq_length = len(event_bases)
    strand = aln_data['mapped_strand']

    kmers_signal = []
    pattern = re.compile(kmer_filter)
    for idx in range(2, seq_length - 3):
        bases = [event_bases[idx + x].decode() for x in [-2, -1, 0, 1, 2]]
        kmer = ''.join(bases)

        if not pattern.search(kmer):
            continue

        if strand == '+':
            pos = aln_data['mapped_start'] + idx + 1  # 1-based coordinate
        else:
            pos = aln_data['mapped_end'] - idx

        encoded_kmer = onehot_encode(kmer)
        encoded_kmer_flat = [val for row in encoded_kmer for val in row]

        mean = [event_mean[idx + x] for x in [-2, -1, 0, 1, 2]]
        std = [event_stdev[idx + x] for x in [-2, -1, 0, 1, 2]]
        length = [event_lengths[idx + x] for x in [-2, -1, 0, 1, 2]]
        chrom = [aln_data['mapped_chrom']]

        feature = chrom + [pos] + [strand] + [read_id] + [idx] + [kmer] + mean + std + length + encoded_kmer_flat
        kmers_signal.append(feature)

    return kmers_signal


def extract_file(file):
    """Wrapper to extract features from a single file"""
    kmer_filter = f"[ACTG][ACTG]{FLAGS.center_base}[ACTG][ACTG]"
    try:
        read_id, aln_data, corr_data = extract_read_signals(file)
        signals = extract_kmer_signals(read_id, aln_data, corr_data, kmer_filter)
    except Exception as e:
        print(str(e))
        return None

    return signals


def iterate_files(file_list):
    """Process files in parallel and collect features into a DataFrame"""
    results = []
    pool = multiprocessing.Pool(processes=int(FLAGS.num_process))
    for file in file_list:
        result = pool.apply_async(extract_file, (file,))
        results.append(result)
    pool.close()

    pbar = tqdm(total=len(file_list), position=0, leave=True)
    df = pd.DataFrame()
    for result in results:
        feature = result.get()
        if feature:
            df = pd.concat([df, pd.DataFrame(feature)], axis=0)
        pbar.update(1)

    pool.join()
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract read-level 5-mer features')
    parser.add_argument('-o', '--output', required=True, help='Feature output directory')
    parser.add_argument('-t', '--num_process', default=1, help='Number of CPUs to use')
    parser.add_argument('-i', '--input', required=True, help='Input fast5 directory')
    parser.add_argument('--group', default='RawGenomeCorrected_000', help='tombo suffix')
    parser.add_argument('--center_base', required=True, help='Central base of 5-mer filter (A/T/C/G)')

    args = parser.parse_args(sys.argv[1:])
    global FLAGS
    FLAGS = args

    # Collect all fast5 file paths
    total_fl = []
    for current_dir, subdirs, files in os.walk(FLAGS.input):
        for filename in files:
            if filename.endswith('.fast5'):
                absolute_path = os.path.abspath(os.path.join(current_dir, filename))
                total_fl.append(absolute_path.rstrip())

    os.makedirs(FLAGS.output + '/tmp', exist_ok=True)

    # Process files in batches of 5000
    for i in range(0, len(total_fl), 5000):
        df = iterate_files(total_fl[i:i + 5000])
        df.to_csv(f"{FLAGS.output}/tmp/{i}.csv", index=False)
