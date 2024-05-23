# python code that takes all the fna or fasta files within a folder and 
# creates a table with the contig number per genome

import os
import sys
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from tqdm import tqdm

# Path to the folder with the fna or fasta files
path = sys.argv[1]

# function that takes the path to the folder and returns a table with the contig number per genome
def contig_count(path):
    # create an empty list to store the results
    results = []
    # loop through all the files in the folder
    for file in tqdm(os.listdir(path), desc="Processing files", unit="file"):
        if file.endswith(".fna") or file.endswith(".fasta"):
            # read the file
            records = list(SeqIO.parse(os.path.join(path, file), "fasta"))
            # get the number of contigs
            contig_number = len(records)
            # get the genome name
            genome_name = file.split(".")[0]
            # append the results to the list
            results.append([genome_name, contig_number])

    # create a pandas dataframe from the results
    df = pd.DataFrame(results, columns=["Genome", "Contig_number"])

    # save the dataframe to a csv file
    df.to_csv("contig_count.csv", index=False)

    return df

# function to plot a histogram of the contig number
def plot_histogram(df):

    plt.hist(df["Contig_number"], bins=50)
    plt.xlabel("Contig number")
    plt.ylabel("Number of genomes")
    plt.title("Histogram of contig number per genome")
    plt.show()
    plt.savefig("contig_histogram.png")

# run the functions
if __name__ == "__main__":
    print("Running contig count\n")
    df = contig_count(path)
    plot_histogram(df)
