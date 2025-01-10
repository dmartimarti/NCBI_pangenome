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
echo "Ancestor_GO,Query_GO" > "$output_file"

# create tmp folder
mkdir -p tmp

# Read GO terms into an array
IFS=$'\n' read -r -d '' -a go_terms < <(tr -d '\r' < "$go_term_file")

# Function to process each GO term
process_go_term() {
    local go_term=$1
    local obo_file=$2
    local owltools_path=$3
    local output_file=$4

    go_term=$(echo "$go_term" | tr -d '\n')
    echo "Processing GO term: $go_term"

    # Run OWLtools and save output in a tmp file
    "$owltools_path" "$obo_file" --ancestors "$go_term" > tmp/$go_term.txt

    # get GO terms from the tmp file, sort, and remove duplicates
    echo "Extracting GO terms from $go_term.txt"
    sed 's/.*\(GO:[0-9]*\).*/\1/' < tmp/$go_term.txt > tmp/${go_term}_filt.txt

    # Append the descendant and query GO term to the output CSV file
    cat tmp/${go_term}_filt.txt | while read descendant; do
        echo "$descendant,$go_term" >> "$output_file"
    done
}

export -f process_go_term

# Run the process_go_term function in parallel using GNU Parallel and pv for progress monitoring
total_go_terms=${#go_terms[@]}
parallel -j 6 process_go_term ::: "${go_terms[@]}" ::: "$obo_file" ::: "$owltools_path" ::: "$output_file" | pv -l -s $total_go_terms > /dev/null

# Remove tmp folder
rm -r tmp

echo "Descendants saved to $output_file"