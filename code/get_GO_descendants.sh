#!/bin/bash

# Usage: ./get_descendants.sh <go_term_file> <obo_file> <output_file>

# Check if the correct number of arguments is provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <go_term_file> <obo_file> <output_file>"
    exit 1
fi

# Assign input arguments to variables
go_term_file=$1
obo_file=$2
output_file=$3

# Path to the owltools executable
owltools_path="/Users/danmarti/Documents/MRC_postdoc/My_projects/NCBI_pangenome/owltools/OWLTools-Runner/target/owltools" 

# Create the output CSV file with headers
echo "Descendant_GO,Query_GO" > "$output_file"

# Iterate through each GO term in the input file
while IFS= read -r go_term; do
    echo "Processing GO term: $go_term"

    # Run OWLtools to get descendants of the GO term, filter for GO terms, extract GO IDs, sort, and remove duplicates
    echo "Running command: $owltools_path $obo_file --descendants $go_term"
    descendants=$("$owltools_path" "$obo_file" --descendants "$go_term" | sed 's/.*\(GO_[0-9]*\).*/\1/' | sort | uniq)

    # Iterate through each descendant GO term
    for descendant in $descendants; do
        # Replace "_" with ":" in the descendant GO term
        descendant=$(echo "$descendant" | sed 's/_/:/g')
        # Replace "_" with ":" in the query GO term
        go_term=$(echo "$go_term" | sed 's/_/:/g')

        # Append the descendant and query GO term to the output CSV file
        echo "$descendant,$go_term" >> "$output_file"
    done
done < "$go_term_file"

echo "Descendants saved to $output_file"