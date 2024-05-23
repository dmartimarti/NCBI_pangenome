# Main bash script for HPC jobs

# root folder for the project
/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli

## Download the genomes using ncbi datasets package
cd /rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/tables

module load anaconda3/personal
source activate ncbi_datasets

datasets download genome accession --inputfile NCBI_ids.txt --filename e_coli_NCBI_genomes.zip --assembly-version latest


### pbs script to download the genomes --------------------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=4:mem=16gb

# Load modules for any applications

module load anaconda3/personal
source activate ncbi_datasets

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/tables

datasets download genome accession --inputfile $WORK/NCBI_ids.txt --filename $WORK/e_coli_NCBI_genomes.zip --assembly-version latest

### --------------------------

unzip e_coli_NCBI_genomes.zip

mv ncbi_dataset/data/*/*.fna ./
rm -r ncbi_dataset
rm README.md

# create 10 folders
for i in {1..9}; do mkdir genomes_$i; done

# move the genomes into the folders

# ls -d -1 $PWD/*.fna > genome_list.txt

# after this step I have moved all genomes into a folder called genomes: /rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes
cd /rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes

# for i in {1..9}; do ls -d -1 $PWD/genomes_$i/*.fna > genome_list_$i.txt; done

ls -d -1 $PWD/genomes_assemblies/*.fna > genome_list.txt

### Bakta annotation --------------------------

mkdir /rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes/bakta_annot/gff3_files

## REMEMBER TO CHANGE VARIABLE NAMES AND LIST LENGTH FOR -J

#PBS -l walltime=3:0:0
#PBS -l select=1:ncpus=24:mem=24gb
#PBS -J 1-9572 

### script 

# Load modules for any applications

module load anaconda3/personal
source activate bakta

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes
INPUT=$WORK/genome_list.txt

GENOME=`head -n $PBS_ARRAY_INDEX $INPUT | tail -1`
GENOME_NAME=$(basename $GENOME)  

bakta --db /rds/general/project/lms-cabreiro-analysis/live/bakta_db/db/ --verbose --output $WORK/bakta_annot/$GENOME_NAME --force --prefix $GENOME_NAME --threads 24 $GENOME
mv $WORK/bakta_annot/$GENOME_NAME/*.gff3 $WORK/bakta_annot/gff3_files/$GENOME_NAME.gff3
rm -r $WORK/bakta_annot/$GENOME_NAME

###
# change every .fasta file into a .fna file


## calculate pairwise mash distances for all genomes

#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=12:mem=24gb

module load anaconda3/personal
source activate panaroo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes

mash sketch -p 12 -o $WORK/mash_sketches $WORK/genomes_*/*.fna
mash dist -p 12 $WORK/mash_sketches.msh $WORK/mash_sketches.msh > $WORK/mash_distances.tab


### filter sequences with differences higher than 0.05
# in the file mash_distances.tab, filter the lines where the third column is higher than 0.05
awk '$3 > 0.05' mash_distances.tab > mash_distances_filtered.tab


## classify phylogroups with EZ_Clermont

#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=24:mem=24gb

module load anaconda3/personal
source activate ezclermont_env

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes

ls $WORK/genomes_*/*.fna | parallel -j24 "ezclermont {} 1>> $WORK/phylogroups_ezclermont.txt  2>> $WORK/results.log"




###
# panaroo -i *.gff3 -o results --clean-mode strict -t 2 --remove-invalid-genes

# python contig_count.py /rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes/genomes_2 && mv contig_count.csv contig_genomes_2.csv
#