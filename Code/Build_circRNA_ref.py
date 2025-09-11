import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def csv_to_fasta(csv_file: str, ref_fasta: str, output_fasta: str):
    """
    Convert a CSV file with genomic coordinates to a FASTA file.
    CSV format (no header): pos_id, chrom, start, end, strand
    """
    # Load reference genome
    seq_records = {record.id: record.seq for record in SeqIO.parse(ref_fasta, "fasta")}

    # Read CSV
    df = pd.read_csv(csv_file, header=None)
    df[[2, 3]] = df[[2, 3]].astype(int)
    df[2] = df[2] - 1  # convert start to 0-based

    fasta_records = []

    for row in df.itertuples(index=False):
        try:
            pos_id = str(row[0])   # unique identifier
            chrom = row[1]
            start = int(row[2])
            end = int(row[3])
            strand = row[4]

            seq = seq_records[chrom][start:end]
            if strand == '-':
                seq = seq.reverse_complement()

            merged_seq = seq + seq
            fasta_records.append(SeqRecord(merged_seq, id=pos_id, description=""))
        except Exception as e:
            print(f"Error processing row {row}: {e}")

    # Write to FASTA
    SeqIO.write(fasta_records, output_fasta, "fasta")
    print(f"FASTA written to {output_fasta}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert CSV to FASTA")
    parser.add_argument("-c", "--csv_file", required=True, help="Path to input CSV file")
    parser.add_argument("-r", "--ref", required=True, help="Path to reference genome FASTA")
    parser.add_argument("-o", "--output_fasta", required=True, help="Path to output FASTA file")

    args = parser.parse_args()
    csv_to_fasta(args.csv_file, args.ref, args.output_fasta)
