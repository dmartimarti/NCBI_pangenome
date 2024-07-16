# python script that reads a multifasta file and splits it into single fasta files

import argparse
import os
from Bio import SeqIO
from tqdm import tqdm

# arguments
parser = argparse.ArgumentParser(description='Split a multifasta file into single fasta files')
parser.add_argument('-i',
                    '--fasta_file',
                    required=True,
                    help='Input multifasta file')

parser.add_argument('-o',
                    '--output_dir',
                    required=True,
                    help='Output directory')

args = parser.parse_args()

# create output directory

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

def split_fasta(fasta_file, output_dir):
    for record in tqdm(SeqIO.parse(fasta_file, "fasta"), total=len(list(SeqIO.parse(fasta_file, "fasta")))):
        output_file = os.path.join(output_dir, record.id + ".fasta")
        with open(output_file, "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")
    return

def main():
    split_fasta(args.fasta_file, args.output_dir)
    return

if __name__ == '__main__':
    main()