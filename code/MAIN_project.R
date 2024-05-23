
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
library(cowplot)
library(viridis)
library(here)
library(openxlsx)
library(patchwork)

theme_set(theme_cowplot(14))


# helper functions -------------------

show_in_excel = function(.data) {
  if (interactive()) {
    tmp = tempfile(fileext = ".csv")
    readr::write_excel_csv(.data, tmp)
    
    fs::file_show(tmp)
  }
  .data
}


# filter metadata -----------------------------------------------------------

# only need to do this once

# metadata = read_xlsx('ecoli_ncbi_metadata.xlsx') 
# 
# metadata = metadata %>% 
#   filter(assembly_level %in% c("Complete Genome",
#                                  "Chromosome",
#                                  "Scaffold")) %>% 
#   filter(checkM_marker_set == 'Escherichia coli')
# 
# 
# 
# 
# write.xlsx(metadata, "ecoli_ncbi_metadata.xlsx",
#            overwrite = TRUE)


## read metadata -----------------------------------------------------------
## START FROM HERE 

metadata = read_xlsx('ecoli_ncbi_metadata.xlsx') 

metadata = metadata %>% 
  # remove the string after the dot in assembly_id
  mutate(assembly_id_simp = str_remove(assembly_id, "\\..*"),
         .before = "assembly_name")


metadata %>% 
  count(assembly_level) %>% 
  print

# unique metadata identifiers 
metadata = metadata %>% 
  distinct(assembly_name, .keep_all = T)
  
metadata %>% 
  filter(checkM_completeness > 90,
         checkM_contamination < 1) 





## plot metadata -----------------------------------------------------------

# plot distribution of contig_n50
metadata %>% 
  ggplot(aes(contig_n50)) +
  geom_histogram(binwidth = 10000) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Contig N50", y = "Number of genomes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis(discrete = T)

ggsave("exploration/contig_n50_distribution.pdf", width = 6, height = 4)

# plot distribution of scaffold_n50
metadata %>% 
  ggplot(aes(scaffold_n50)) +
  geom_histogram(binwidth = 10000) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Scaffold N50", y = "Number of genomes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("exploration/scaffold_n50_distribution.pdf", width = 6, height = 4)

metadata %>% 
  filter(contig_n50 > 4000000)


# plot scaffold_n50 vs contig_n50
metadata %>% 
  ggplot(aes(contig_n50, scaffold_n50)) +
  geom_point(alpha = 0.6) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Contig N50", y = "Scaffold N50") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1)) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(method = "pearson", size = 4)

ggsave("exploration/contig_vs_scaffold_n50.pdf", width = 6, height = 4)

# plot scaffold distribution by assembly level with violin plot
metadata %>% 
  ggplot(aes(assembly_level, scaffold_n50, fill = assembly_level)) +
  geom_violin(show.legend = F) +
  # plot mean and sd as pointrange
  stat_summary(fun.data = mean_sdl, geom = "pointrange", 
               fun.args = list(mult = 1), color = "black",
               show.legend = F) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Assembly level", y = "Scaffold N50") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis(discrete = T, begin = 0.2, end = 0.8) +
  theme_cowplot(14)

ggsave("exploration/scaffold_n50_by_assembly_level.pdf", width = 6, height = 4)

# show scaffold distribution of assembly_level == scaffold
p1 = metadata %>% 
  filter(assembly_level == 'Scaffold') %>% 
  ggplot(aes(scaffold_n50)) +
  geom_histogram(binwidth = 10000) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Scaffold N50", y = "Number of genomes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_cowplot(14)

p2 = metadata %>% 
  filter(assembly_level == 'Scaffold') %>%
  filter(scaffold_n50 < 500000) %>% 
  ggplot(aes(scaffold_n50)) +
  geom_histogram(binwidth = 10000) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Scaffold N50", y = "Number of genomes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_cowplot(14)


# show them with patchwork
p1 + p2

ggsave("exploration/scaffold_n50_distribution_scaffold.pdf", width = 8, height = 5)

metadata_filt = metadata %>% filter(scaffold_n50 > 100000)
dims = dim(metadata_filt)
dims[1]



## checkM ------------------------------------------------------------------



# empty values from checkM 

metadata_filt %>% 
  mutate(checkM_completeness_class = 
           case_when(is.na(checkM_completeness) ~ 'NA',
                     checkM_completeness < 95 ~ "< 95%",
                     checkM_completeness >= 95 ~ ">= 95%")) %>% 
  count(checkM_completeness_class) %>% print

metadata_filt %>% 
  mutate(checkM_contamination_class = 
           case_when(is.na(checkM_contamination) ~ 'NA',
                     checkM_contamination > 1 ~ "> 1%",
                     checkM_contamination <= 1 ~ "<= 1%")) %>% 
  count(checkM_contamination_class) %>% print


dim_filt = metadata_filt %>%
  filter(checkM_completeness > 95,
         checkM_contamination < 1) %>% dim


# version of the metadata after filtering for checkM
metadata_filt = metadata_filt %>% 
  filter(checkM_completeness > 95,
         checkM_contamination < 1)



## filtering extreme cases of genome size ------------------------

metadata_filt %>% 
  # plot distribution of seq_length 
  ggplot(aes(seq_length)) +
  geom_histogram(binwidth = 10000) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Sequence length", y = "Number of genomes")
# filter genomes with seq_length > 7000000

ggsave("exploration/seq_length_distribution.pdf", width = 6, height = 4)

metadata_filt %>% 
  # plot distribution of genes 
  ggplot(aes(genes)) +
  geom_histogram() +
  # scale_x_continuous(labels = scales::comma) +
  labs(x = "Genes", y = "Number of genomes")
# filter genomes with genes > 7000 and genes < 2000

ggsave("exploration/genes_distribution.pdf", width = 6, height = 4)


# filter genomes with seq_length > 7000000
metadata_filt = metadata_filt %>% 
  filter(seq_length < 7000000)

  
# write_csv(metadata_filt, "ecoli_ncbi_metadata_filt.csv")
write_csv(metadata_filt, "tables/ecoli_ncbi_metadata_filt.csv")


## get assembly names
metadata_filt %>% 
  select(assembly_id) %>% 
  write_csv("tables/assembly_id.csv")
  


#  read contig files ------------------------------------------------------
contigs = fs::dir_ls("data/genome_contigs") %>% # list all files 
  map_df(read_csv, .id = "file") %>% # read all files and bind them together
  rename(assembly_id_simp = Genome) %>% 
  select(-file)

contigs %>% write_csv("tables/genome_contigs.csv")

# read all csv files within data/genome_contigs and bind them together
metadata_filt =  contigs %>% 
  left_join(metadata_filt)

# plot histogram of number of contigs

metadata_filt %>% 
  filter(Contig_number < 500) %>% 
  ggplot(aes(Contig_number)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Number of contigs", y = "Number of genomes") +
  theme_cowplot(14)

ggsave("exploration/contig_number_distribution.pdf", width = 6, height = 4)

metadata_filt = metadata_filt %>% 
  filter(Contig_number < 300) 

# update the metadata_filt files
# write_csv(metadata_filt, "ecoli_ncbi_metadata_filt.csv")
write_csv(metadata_filt, "tables/ecoli_ncbi_metadata_filt.csv")

# create a list of genomes with high contigs to remove from analysis
contigs %>% 
  filter(Contig_number > 300)  %>% 
  write_csv("tables/REMOVALS_high_contigs.csv")



# mash distances ----------------------------------------------------------

# mash distances have been calculated pairwise between all genomes
# after that I have filtered those distances that were > 0.05

# get the mash distances

mash_dist = read_tsv("tables/mash_distances_filtered.tsv",
         col_names = F) %>% 
  rename(query_genome = `X1`,
         target_genome = `X2`,
         mash_distance = `X3`)


# plot distribution of mash distances
mash_dist %>% 
  ggplot(aes(mash_distance)) +
  geom_histogram(binwidth = 0.001) +
  labs(x = "Mash distance", y = "Number of pairs") +
  theme_cowplot(14)

mash_dist %>% 
  count(query_genome) %>% 
  arrange(desc(n))



## phylogroups -----------------------------------------------------------

phylogroups = read_tsv("tables/phylogroups_ezclermont.tsv",
         col_names = F) %>% 
  rename(genome_id = `X1`,
         phylogroup = `X2`) %>% 
  mutate(assembly_id = str_sub(genome_id, 1, 15))

# update the metadata with the phylogroup information
metadata_filt = metadata_filt %>% 
  left_join(phylogroups) 


# plot a histogram of phylogroups
metadata_filt %>% 
  ggplot(aes(phylogroup, fill = phylogroup)) +
  geom_bar() +
  labs(x = "Phylogroup", y = "Number of genomes") +
  theme_cowplot(14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# FINAL VERSION OF THE METADATA -------------------------------------------

removals_contigs = read_csv("tables/REMOVALS_high_contigs.csv") %>% 
  pull(assembly_id_simp)


removals_phylogroup = metadata_filt %>%
  filter(phylogroup %in% c("cryptic", 
                           "EC_control_fail",
                           "U/cryptic",
                           NA)) %>% 
  pull(assembly_id_simp)


c(removals_contigs, removals_phylogroup) %>% 
  unique() %>% 
  as_tibble() %>% 
  # save a file with the removals 
  write_delim("tables/REMOVALS_final.txt", delim = "\n")


## add the lab genomes ---------------------


lab_metadata =
   read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx")

lab_str = read_csv("tables/lab_strains.txt", 
         col_names = FALSE) %>% 
  rename(fasta = `X1`)

lab_metadata = lab_metadata %>% 
  filter(fasta %in% lab_str$fasta) %>% 
  distinct(fasta, .keep_all = T) %>% 
  select(assembly_ID = ID, 
         assembly_name = Strainname, 
         phylogroup = phylogroup) %>% 
  mutate(assembly_name = case_when(is.na(assembly_name) ~ assembly_ID,
                                   TRUE ~ assembly_name)) 






## filter the metadata and save it
metadata_final = metadata_filt %>% 
  bind_rows(lab_metadata) %>%
  filter(!assembly_id_simp %in% c(removals_contigs, removals_phylogroup))

metadata_final %>% 
  write_csv("tables/MAIN_ecoli_NCBI_metadata.csv")




# Including strains at contig level ---------------------------------------

# I'm thinking that I may have the ability to include strains at the contig level
# I will filter for contigs and then filter for the same criteria as before

ecoli_contigs = read_excel("ecoli_ncbi_metadata_complete.xlsx") %>% 
  filter(assembly_level == "Contig") %>% 
  distinct(assembly_name, .keep_all = T) %>% 
  filter(checkM_completeness > 95,
         checkM_contamination < 1) %>% 
  filter(seq_length < 7000000) 

ecoli_contigs %>%
  select(assembly_id) %>%
  write_delim("tables/contigs/contigs_id.txt", delim = "\n")




# Phylogroup distribution -------------------------------------

phylo_colors = c("A" = "#1b9e77",
                 "B1" = "#d95f02",
                 "B2" = "#7570b3",
                 "D" = "#e7298a",
                 "E" = "#66a61e",
                 "F" = "#e6ab02",
                 "G" = "#a6761d",
                 "C" = "#3A254D",
                 "U" = "#666666")

metadata_final %>% 
  ggplot(aes(phylogroup, fill = phylogroup)) +
  geom_bar(show.legend = F) +
  labs(x = "Phylogroup", y = "Number of genomes") +
  scale_fill_manual(values = phylo_colors) +
  theme_cowplot(14) 

ggsave("exploration/phylogroup_distribution.pdf", width = 6, height = 4)





