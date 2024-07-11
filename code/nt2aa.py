from Bio import SeqIO
from Bio.Seq import Seq
import argparse
from tqdm import tqdm

# Create argument parser
parser = argparse.ArgumentParser(description='Translate nucleotide sequences to amino acid sequences and remove sequences with early stop codons.')
parser.add_argument('-i', '--input', type=str, required=True, help='Input FASTA file path')
parser.add_argument('-o', '--output', type=str, required=True, help='Output FASTA file path')
args = parser.parse_args()

# Open input and output files
with open(args.input, 'r') as input_file, open(args.output, 'w') as output_file:
    # Initialize counter for discarded sequences
    discarded_count = 0
    # Loop through input sequences with tqdm 
    for record in SeqIO.parse(input_file, 'fasta'):
        # Translate sequence to amino acids
        aa_seq = record.seq.translate(table=11, to_stop=True)
        # Check if sequence has early stop codons
        if '*' not in aa_seq:
            # Write sequence to output file
            output_file.write(f'>{record.id}\n{aa_seq}\n')
        else:
            # Increment counter for discarded sequences
            discarded_count += 1

# Print count of discarded sequences
print(f'{discarded_count} sequences were discarded due to early stop codons.')
