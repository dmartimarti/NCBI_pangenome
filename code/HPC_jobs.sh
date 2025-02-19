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

### --------------------------
# compress the bakta folder into a gz file with pigz

#PBS -l walltime=6:0:0
#PBS -l select=1:ncpus=24:mem=24gb

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes/bakta_annot

# compress the folder into a gz file with 20 threads and keep the original folder
tar --use-compress-program="pigz -k -p 20" -cf $WORK/gff3_files.tar.gz $WORK/gff3_files



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


### run PPanGGOLiN

# make a table with the genome names and their path
find $PWD -type f -name "*.fna.gff3" | awk -v OFS="\t" -F/ '{ print substr($NF, 1, length($NF)-9), $0 }' > genomes.gbff.txt


#PBS -l walltime=48:0:0
#PBS -l select=1:ncpus=32:mem=150gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes/bakta_annot/gff3_files
OUT=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli

ppanggolin all --anno $WORK/genomes.gbff.txt --output $OUT/ppanggolin_results -c 32 --verbose 1 -f


###
# panaroo -i *.gff3 -o results --clean-mode strict -t 2 --remove-invalid-genes

# python contig_count.py /rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes/genomes_2 && mv contig_count.csv contig_genomes_2.csv
#
find $PWD -type f -name "*.fna.gff3" | awk -v OFS="\t" -F/ '{ print substr($NF, 1, length($NF)-9), $0 }' > genomes.gbff.txt
ppanggolin all --anno genomes.gbff.txt --output ppanggolin_results -c 2 --verbose 2 -f

# IMPORTANT, READ HERE
# IMPORTANT, READ HERE
# IMPORTANT, READ HERE
# IMPORTANT, READ HERE
# IMPORTANT, READ HERE
# IMPORTANT, READ HERE
# IMPORTANT, READ HERE
# IMPORTANT, READ HERE
# IMPORTANT, READ HERE
# IMPORTANT, READ HERE

# some of the gff files have non-ASCII characters that are causing ppanggolin to crash
# for example:
LC_ALL=C grep -n -P [$'\x80'-$'\xFF'] AUS_30.fna.gff3

# this will show that the file AUS_30.fna.gff3 has non-ASCII characters

# 1776:contig_13	Prodigal	CDS	25409	25765	.	-	0	ID=FBPPNK_08630;Name=putative DNA-binding protein with ???double-wing??? 
# structural motif%2C MmcQ/YjbR family;locus_tag=FBPPNK_08630;product=putative DNA-binding protein with ???double-wing??? structural motif%2C 
# MmcQ/YjbR family;Dbxref=COG:COG2315,COG:K,RefSeq:WP_000153726.1,SO:0001217,UniParc:UPI0002185663,UniRef:UniRef100_A0A066SY50,UniRef:UniRef50_P0AF51,UniRef:UniRef90_A0A234IDG1;
# gene=mmcQ

# gff3 files can be fixed by removing the non-ASCII characters with this command:
for file in *.gff3; do
  LC_ALL=C tr -cd '\0-\177' < "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
done





### ppanggolin with annotation ------------------------

find $PWD -type f -name "*.fna" | awk -v OFS="\t" -F/ '{ print substr($NF, 1, length($NF)-4), $0 }' > genomes.gbff.txt


#PBS -l walltime=72:0:0
#PBS -l select=1:ncpus=32:mem=196gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes/genomes_assemblies
OUT=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli

ppanggolin all --fasta $WORK/genomes.gbff.txt --output $WORK/ppanggolin_annot_results -c 32 --verbose 2 -f



### panaroo script --------------------------

#PBS -l walltime=72:0:0
#PBS -l select=1:ncpus=32:mem=196gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes/bakta_annot/gff3_files
OUT=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes

panaroo -i $WORK/*.gff3 -o $OUT/panaroo_results --clean-mode strict -t 32 --remove-invalid-genes -a core




### PPanGGOLiN graph --------------------------

#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=24:mem=296gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/ppanggolin_results

ppanggolin graph -p $WORK/pangenome.h5




### PPanGGOLiN partition --------------------------

#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=24:mem=296gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/ppanggolin_results

ppanggolin partition -p $WORK/pangenome.h5 -c 24 --draw_ICL -o $WORK/partition_results --verbose 2


### PPanGGOLiN partition --------------------------

#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=12:mem=296gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/ppanggolin_results

ppanggolin rgp -p $WORK/pangenome.h5 -c 12 -o $WORK/rpg_results --verbose 1

ppanggolin spot -p $WORK/pangenome.h5 -c 12 --spot_graph --graph_formats gexf -o $WORK/spot_results --verbose 1


### PPanGGOLiN module --------------------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=12:mem=146gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/ppanggolin_results

ppanggolin module -p $WORK/pangenome.h5 -c 12 --verbose 1






### PPanGGOLiN draw --------------------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=12:mem=96gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/ppanggolin_results

ppanggolin draw -p $WORK/pangenome.h5 --tile_plot --draw_spots --ucurve -o $WORK/draw_results --verbose 1




### PPanGGOLiN stats --------------------------

#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=16:mem=96gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/ppanggolin_results

ppanggolin write_pangenome -p $WORK/pangenome.h5 -c 16 -o $WORK/stats_results --verbose 1 --gexf --json --csv --Rtab  --stats \
    --partitions --families_tsv --regions --spots --modules --spot_modules


### PPanGGOLiN MSA --------------------------

#PBS -l walltime=72:0:0
#PBS -l select=1:ncpus=32:mem=596gb

module load anaconda3/personal
source activate ppanggo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/ppanggolin_results

ppanggolin msa -p $WORK/pangenome.h5 -c 32 -o $WORK/msa --verbose 1 --partition persistent


### panaroo script --------------------------

#PBS -l walltime=72:0:0
#PBS -l select=1:ncpus=64:mem=500gb

module load anaconda3/personal
source activate panaroo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/genomes/bakta_annot/gff3_files
OUT=/rds/general/project/lms-cabreiro-analysis/ephemeral/NCBI_ecoli/panaroo

panaroo -i $WORK/gff3_batch_1/*.gff3 -o $OUT/results_batch_1 --clean-mode strict -t 64 --remove-invalid-genes --merge_paralogs


### Panaroo merge for all results

#PBS -l walltime=72:0:0
#PBS -l select=1:ncpus=64:mem=360gb

# Load modules for any applications

module load anaconda3/personal
source activate panaroo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results

panaroo-merge -d  $WORK/results_batch_1 $WORK/results_batch_2 $WORK/results_batch_3 $WORK/results_batch_4 $WORK/results_batch_5  \
 	-o $WORK/panaroo_final -t 64 




### Mafft alignment of core genomes --------------------------

#PBS -l walltime=48:0:0
#PBS -l select=1:ncpus=32:mem=64gb

# Load modules for any applications

module load anaconda3/personal
source activate panaroo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_no_paralogs/core_genome_genes

/rds/general/user/dmarti14/home/anaconda3/envs/panaroo/bin/mafft --thread 2 --auto --clustalout --reorder $WORK/core_91.fasta > $WORK/core_91.aln


###### INTERPRO ------------

#PBS -l walltime=48:0:0
#PBS -l select=1:ncpus=48:mem=64gb

# Load modules for any applications

module load anaconda3/personal
module load java/oracle-jdk-11.0.10
source activate interpro


WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_paralogs

/rds/general/user/dmarti14/home/interproscan/interproscan.sh -i $WORK/pan_genome_reference_aa_paralogs.fa -cpu 48 -dp -d $WORK -goterms -f tsv




# for paralogs --------------------------
#PBS -l walltime=48:0:0
#PBS -l select=1:ncpus=48:mem=64gb

# Load modules for any applications

module load anaconda3/personal
module load java/oracle-jdk-11.0.10
source activate interpro


WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_paralogs

/rds/general/user/dmarti14/home/interproscan/interproscan.sh -i $WORK/pan_genome_reference_aa.fa -cpu 48 -dp -d $WORK -goterms -f tsv




###### Eggnog ------------

download_eggnog_data.py -F -P -M -H --data_dir /rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/eggnog_db


#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=40:mem=80gb

# Load modules for any applications
module load anaconda3/personal
source activate eggnog

DB=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/eggnog_db
WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_paralogs

emapper.py -i $WORK/pan_genome_reference.fa -o $WORK/complete_PG_eggnogmap --data_dir $DB --override --itype CDS --cpu 40 -m novel_fams -m mmseqs


# for paralogs --------------------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=40:mem=80gb

# Load modules for any applications
module load anaconda3/personal
source activate eggnog

DB=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/eggnog_db
WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_paralogs

emapper.py -i $WORK/pan_genome_reference_paralogs.fa -o $WORK/complete_PG_eggnogmap_paralogs \
--data_dir $DB --override --itype CDS --cpu 40 -m novel_fams -m mmseqs



### PROGENOMES DB --------------------

#!/bin/bash
#PBS -N diamond-db-build
#PBS -l select=1:ncpus=20:mem=200gb
#PBS -l walltime=24:00:00

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_no_paralogs/progenomes

echo "Unziping representative proteins"

bunzip2 $WORK/progenomes3.proteins.representatives.fasta.bz2

# Load DIAMOND
module load anaconda3/personal
source activate progenomes

# Make a local copy of the input file in the node's scratch space
# cp $WORK/proteins.fasta $TMPDIR/progenomes3.proteins.representatives.fasta

# Make the DIAMOND database
echo "### RUNNING DIAMOND MAKEDB"
diamond makedb --in $WORK/progenomes3.proteins.representatives.fasta -d  $WORK/progenome3.proteins --threads 20


### PROGENOMES BLASTX ------------------------------

#!/bin/bash

#PBS -N diamond-align
#PBS -l select=1:ncpus=32:mem=200gb
#PBS -l walltime=24:00:00

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_paralogs
DB=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_no_paralogs/progenomes

# Load DIAMOND
module load anaconda3/personal
source activate progenomes

# Run DIAMOND to align the query nucleotide sequences against the protein database
diamond blastx --query $WORK/pan_genome_reference.fa --db $DB/progenome3.proteins.dmnd --out $WORK/progenomes/reference_PG_progenomes3.txt --threads 32 --outfmt 6




### Proteinfer -----------------------------

#PBS -l walltime=48:00:00,select=1:ncpus=16:mem=64G:ngpus=1:gpu_type=RTX6000

module load anaconda3/personal
module load cuda/8.0
module load cudnn/6.5

source activate proteinfer

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_paralogs/aa_files
MODELS=$HOME/pangenome_study/progenomes/gene_seqs/cached_models

pytho $HOME/proteinfer/proteinfer.py -i $WORK/pan_1.fa -o $WORK/pan_1.tsv --num_ensemble_elements 5 --reporting_threshold 0.5 \
        --model_cache_path $MODELS





# --------------------------
# panaroo extract genes 


#!/bin/bash
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -l walltime=24:00:00

module load anaconda3/personal
source activate panaroo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_no_paralogs




panaroo-extract-gene -q accC --pa $WORK/gene_presence_absence.csv --gene $WORK/gene_data.csv -o $WORK/core_genome_genes --dna --idtype 'isolate' 

# --------------------------


#PBS -l walltime=2:0:0
#PBS -l select=1:ncpus=24:mem=24gb
#PBS -J 1-150 

# Load modules for any applications

module load anaconda3/personal
source activate panaroo

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_no_paralogs/core_genome_genes/extracted_gene_sequences_deduplicated
INPUT=$WORK/gene_list.txt

GENOME=`head -n $PBS_ARRAY_INDEX $INPUT | tail -1`
GENOME_NAME=$(basename $GENOME)  

/rds/general/user/dmarti14/home/anaconda3/envs/panaroo/bin/mafft --thread 24 --auto --clustalout --reorder $WORK/$GENOME > $WORK/aln/$GENOME_NAME.aln






######## PHYLOGENETIC TREE ----------------------------------------

#PBS -l walltime=72:0:0
#PBS -l select=1:ncpus=32:mem=240gb

# Load modules for any applications

module load anaconda3/personal
source activate iqtree

WORK=/rds/general/project/lms-cabreiro-analysis/live/NCBI_ecoli/panaroo_results/panaroo_no_paralogs/core_genome_genes

mkdir $WORK/phy_tree

cp $WORK/core_150.aln $TMPDIR

iqtree -s core_150.aln -fast -m GTR+I+G -T AUTO --threads-max 32 32 -v -pre core_tree

cp core_tree* $WORK/phy_tree




#### transform clusters from CD-HIT to  a tsv file
awk '/^>Cluster/ {cluster=$2} /^>/ {next} {print cluster, $1, $3}' pan_genome_reference_aa_FINAL_unique.fa.clstr | tr ',' '\t' > PCA_embeddings_clusters.tsv




### Bio-embeddings test on the HPC -----------------------------

# EXAMPLE OF CONFIG FILE
global:
  sequences_file: /rds/general/user/dmarti14/home/bio_embeddings/PG_data/pan_genome_reference_aa_FINAL.fa
  prefix: PG_embeddings_UMAP

t5_embeddings:
  type: embed
  protocol: prottrans_t5_xl_u50
  device: cuda
  half_precision_model: True
  half_precision: True
  reduce: True
  discard_per_amino_acid_embeddings: True

umap_projections:
  type: project
  protocol: umap
  depends_on: t5_embeddings
  n_components: 2
  n_jobs: 12




#PBS -l select=1:ncpus=12:mem=64gb:ngpus=1:gpu_type=RTX6000 -lwalltime=1:00:00

module load anaconda3/personal
module load cuda/10.2  
module load cudnn/8.2.4 

source activate bioembeddings

WORK=$HOME/bio_embeddings

bio_embeddings $WORK/PG_data/config_PG_embedT5_UMAP.yml

# remember to change the output folder name in the config file
cp -r $TMPDIR/PG_embeddings_UMAP $WORK