import os
from Bio import SeqIO
import sys

"""
This script reads all the fasta files in the input folder and removes duplicate sequences.
The unique sequences are saved in the output folder with the same filename.

The script is run after the genes have been extracted from the genomes in panaroo. 

It prints the final number of sequences, that must be: 9558
"""

def remove_duplicates_and_save(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta"):
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename)

            # Read sequences and remove duplicates
            unique_sequences = {}
            for record in SeqIO.parse(input_path, "fasta"):
                if record.id not in unique_sequences:
                    unique_sequences[record.id] = record

            # Print the total number of unique sequences
            print(f"File: {filename}, Total Unique Sequences: {len(unique_sequences)}")

            # Write the unique sequences to the output file
            with open(output_path, "w") as output_handle:
                SeqIO.write(unique_sequences.values(), output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python remove_duplicates.py <input_folder> <output_folder>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    remove_duplicates_and_save(input_folder, output_folder)
