import os
import pandas as pd
import argparse
import sys
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import mappy as mp

# 开始处理每个 CSV 文件
def process_directory(directory_path, seq_records, reads_file, output_file, merged_fasta_name):
    csv_files = [f for f in os.listdir(directory_path) if f.endswith('.csv')]
    with tqdm(total=len(csv_files), desc="Total CSV files") as total_pbar:
        for filename in csv_files:
            file_path = os.path.join(directory_path, filename)
            print(f"Processing file: {file_path}", flush=True)

            try:
                df = pd.read_csv(file_path, header=None)
                df[[2, 3]] = df[[2, 3]].astype(int)
                df[2] = df[2] - 1

                fasta_records = []
                seq_strand_map = {}

                for row in df.itertuples(index=False):
                    try:
                        pos_id = row[0]
                        chrom = row[1]
                        start = int(row[2])
                        end = int(row[3])
                        strand = row[4]

                        seq = seq_records[chrom][start:end]
                        if strand == '-':
                            seq = seq.reverse_complement()
                        merged_seq = seq + seq
                        fasta_records.append(SeqRecord(merged_seq, id=pos_id, description=""))
                        seq_strand_map[pos_id] = strand
                    except Exception as e:
                        print(f"Error processing row {row}: {e}")

                fasta_name = merged_fasta_name
                SeqIO.write(fasta_records, fasta_name, "fasta")

                try:
                    aligner = mp.Aligner(fasta_name)
                    if not aligner:
                        raise Exception("Failed to load reference.")

                    print("Aligner initialized", flush=True)
                    results = []
                    for name, seq, qual in mp.fastx_read(reads_file):
                        for hit in aligner.map(seq, cs=True, MD=True):
                            if hit.is_primary and hit.strand == 1:
                                try:
                                    ctg, rng = hit.ctg.split(':')
                                    start, end = map(int, rng.split('-'))
                                    real_start = hit.r_st + start
                                    real_end = hit.r_en + start
                                    if real_start < end - 10 and real_end > end + 10:
                                        original_strand = seq_strand_map.get(hit.ctg, '+')
                                        results.append((name, hit.ctg, original_strand))
                                except Exception as e:
                                    print(f"Error in alignment parsing: {e}")

                    print(f"Found {len(results)} results in {filename}", flush=True)
                    with open(output_file, "a+") as f:
                        for r in results:
                            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")

                except Exception as e:
                    print(f"Alignment error in file {filename}: {e}")

            except Exception as e:
                print(f"Error processing file {filename}: {e}")
            total_pbar.update(1)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find circRNA')
    parser.add_argument('-i', '--input', required=True, help='input fastq file')
    parser.add_argument('-d', '--directory', required=True, help='directory path with CSV files')
    parser.add_argument('-r', '--ref', required=True, help='reference fasta file')
    parser.add_argument('-o', '--output', required=True, help='output file')
    parser.add_argument('-t', '--threads', default=1, type=int, help='number of threads')
    parser.add_argument('-m', '--merged_fasta', default='merged_sequences.fasta', help='output fasta file for merged sequences')

    args = parser.parse_args()

    # 加载参考序列
    seq_records = {record.id: record.seq for record in SeqIO.parse(args.ref, "fasta")}

    process_directory(
        directory_path=args.directory,
        seq_records=seq_records,
        reads_file=args.input,
        output_file=args.output,
        merged_fasta_name=args.merged_fasta
    )

