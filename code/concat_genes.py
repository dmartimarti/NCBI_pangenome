import os
import argparse
from Bio import SeqIO
from tqdm import tqdm

# usage: python code/concat_genes.py data/genes data/concatenated_genes.fasta

def concatenate_genes_from_folder(folder_path, output_file):
    # Dictionary to store concatenated sequences for each genome
    concatenated_sequences = {}
    
    # List all fasta files in the folder
    fasta_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.aln')]
    print(f"Found {len(fasta_files)} fasta files in the folder.")
    for fasta_file in tqdm(fasta_files, total=len(fasta_files)):
        with open(fasta_file, 'r') as handle:
            for record in tqdm(SeqIO.parse(handle, "fasta")):
                genome_id = record.id
                sequence = str(record.seq)
                
                if genome_id not in concatenated_sequences:
                    concatenated_sequences[genome_id] = sequence
                else:
                    concatenated_sequences[genome_id] += sequence
    
    # Write concatenated sequences to output file
    with open(output_file, 'w') as output_handle:
        for genome_id, concatenated_seq in concatenated_sequences.items():
            output_handle.write(f">{genome_id}\n")
            output_handle.write(f"{concatenated_seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Concatenate gene sequences from multiple fasta files into a single file.")
    parser.add_argument('folder_path', type=str, help="Path to the folder containing fasta files.")
    parser.add_argument('output_file', type=str, help="Output file name for concatenated sequences.")
    
    args = parser.parse_args()
    print(args)
    concatenate_genes_from_folder(args.folder_path, args.output_file)

if __name__ == "__main__":
    main()