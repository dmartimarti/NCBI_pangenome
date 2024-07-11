import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool

def process_data(args):
    i, n, df, genomes_subset = args
    genes_subset = df[genomes_subset]
    genes_subset = genes_subset[genes_subset.sum(axis=1) > 0]
    ngenes = len(genes_subset)
    data = [ngenes, n, i]
    return data

def main():
    parser = argparse.ArgumentParser(description="Permutation and data processing script")
    parser.add_argument("-i", "--input", required=True, help="Input file path")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument("-perm", "--permutations", type=int, required=True, help="Number of permutations")
    parser.add_argument("-points", "--points", type=int, required=True, help="Number of points to calculate")
    parser.add_argument("-p", "--processes", type=int, default=4, help="Number of processes to use")
    args = parser.parse_args()

    print(f'Reading the pangenome data into a DataFrame')
    df = pd.concat([chunk for chunk in tqdm(pd.read_csv(args.input, chunksize=10000, sep='\t', index_col=0), desc='Loading data')])
    # df = pd.read_csv(args.input, sep='\t', index_col=0)
    # df = dd.read_csv(args.input, sep="\t", sample = 2000000).set_index("Gene").compute()
    
    nperm = args.permutations

    total_genomes = len(df.columns)

    print("Removing genes present in 0.99 or more of the genomes")
    df = df[df.sum(axis=1) < 0.99 * total_genomes]

    npoints = args.points
    genomes = np.round(np.linspace(1, total_genomes, npoints)).astype(int)

    acc_df = pd.DataFrame(columns=["n_genes", "n_genomes", "permutation"])

    # create a list of arguments for each permutation and number of genomes
    args_list = []
    for i in range(nperm):
        for n in genomes:
            genomes_subset = np.random.choice(df.columns, n, replace=False)
            args_list.append((i, n, df, genomes_subset))

    # create a Pool object with the number of processes
    print(f"Running with {args.processes} processes")
    with Pool(processes=args.processes) as pool:
        # apply the process_data function to each argument in parallel
        results = pool.map(process_data, args_list)

    # create a dataframe from the results
    acc_df = pd.DataFrame(results, columns=["n_genes", "n_genomes", "permutation"])
    
    acc_df.to_csv(args.output, index=False)

    print("Done!")

if __name__ == "__main__":
    main()