
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
library(cowplot)
library(viridis)
library(here)
library(openxlsx)
library(patchwork)
library(ggExtra)
library(DescTools)
library(rstatix)
library(openxlsx)
library(glue)
library(styler)

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
  filter(assembly_level == "Scaffold") %>%
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

# metadata_final %>%
#   write_csv("tables/MAIN_ecoli_NCBI_metadata.csv")

# metadata_final = read_csv("tables/MAIN_ecoli_NCBI_metadata.csv")



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





# PANGENOME ---------------------------------------------------------------


# IMPORTANT VARS ----------------------------------------------------------

# genome_names: tibble with genome names
# metadata_final: metadata with phylogroups
# phylogroups: phylogroups of the genomes
# phylo_colors: colors for the phylogroups
# gene_classes_phylogroup: gene classes per phylogroup
# gene_presence: gene groups (core, shell, cloud)
# refseqs_n: number of reference sequences

# Gene P/A ----------------------------------------------------------------

# read the gene presence absence matrix

gene_pa = read_delim("data/panaroo_results/gene_presence_absence.Rtab", 
           delim = "\t", escape_double = FALSE, 
           trim_ws = TRUE)

names(gene_pa)[2:length(names(gene_pa))] %>% 
  as_tibble() %>%
  mutate(value = str_replace(value, ".fna", ""),
         # remove anything after the first dot when the string starts with GCA
         value = str_remove(value, "\\..*")
         ) %>%
  # save a txt
  write_delim("tables/genome_names.txt", delim = "\n")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ## 
## IMPORTANT: I HAVE FIXED THE NAMES, LOAD THEM FROM THIS FILE ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ## 

genome_names = read_csv("tables/genome_names.txt")

# transpose the data frame
gene_pa_t = gene_pa %>% 
  column_to_rownames("Gene") %>%
  t() %>%
  as_tibble() %>% 
  mutate(genome = genome_names$value,
         .before = "yeeA~~~nadB")

# change the column names to genome names
names(gene_pa) = c("Gene", gene_pa_t$genome)


gene_pa_t %>%
  write_csv("tables/gene_pa_matrix.csv")

# save in Rtab format

gene_pa %>% 
  write_delim("tables/gene_pa_matrix.Rtab", delim = "\t")


## cummulative gene presence absence

gene_presence = gene_pa %>% 
  select(-Gene) %>% 
  as.matrix() %>% 
  rowSums()


gene_presence  = tibble(Gene = gene_pa$Gene,
       gene_presence = gene_presence) %>% 
  mutate(gene_prop = gene_presence / (ncol(gene_pa) - 1))

gene_presence %>% 
  write_csv("tables/gene_presence_proportion.csv")

gene_presence %>% 
  ggplot(aes(gene_prop)) +
  geom_histogram(binwidth = 0.01, fill = "darkblue") +
  labs(
    x = "Proportion of genomes with gene",
    y = "Number of genes"
  ) 

ggsave("exploration/gene_presence_proportion.pdf", width = 7, height = 5)


## PCA of P/A matrix -----------------------------------------------------

# the PCA has been calculated in Python and saved 

# load the PCA results

pca_data = read_delim("tables/PCA_gene_pa_matrix.tsv", 
           delim = "\t", escape_double = FALSE, 
           trim_ws = TRUE) %>% 
  rename(assembly_id_simp = `...1`) %>% 
  # from the assembly_id_simp that starts with "NT", remove all after an underscore "_" if there is one
  mutate(assembly_id_simp = case_when(str_detect(assembly_id_simp, "NT12") ~ str_remove(assembly_id_simp, "_.*"),
                                      TRUE ~ assembly_id_simp)) %>%
  left_join(metadata_final)

metadata_final %>% filter(is.na(phylogroup))


pca_data %>% 
  drop_na(assembly_id_simp) %>% 
  left_join(phylogroups) %>%
  ggplot(aes(PC1, PC2, color = phylogroup)) +
  geom_point() 


pca_data %>% 
  ggplot(aes(PC1, PC2, color = phylogroup,
             fill = phylogroup)) +
  stat_ellipse(type = 't',
               geom  = 'polygon',
               alpha = 0.2) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = phylo_colors) +
  scale_fill_manual(values = phylo_colors) +
  ggrepel::geom_text_repel(
    data = phylo_centroids,
    aes(label = phylogroup,
        x = PC1_centroid,
        y = PC2_centroid),
    color = 'black',
    box.padding = 0.35
  ) +
  labs(
    x = "PC1 [10.87%]",
    y = "PC2 [4.86%]"
  )


ggsave('exploration/PCA_genePA_phylogroups.pdf',
       height = 7, width = 9)




# PG summary stats piechart -----------------------------------------------


PG_summary = read_delim("data/panaroo_results/summary_statistics.txt", 
                        delim = "\t", escape_double = FALSE, 
                        col_names = FALSE, trim_ws = TRUE) %>% 
  rename(Category = X1, 
         Description = X2,
         n_genes = X3) 


# plot piechart

PG_summary %>% 
  filter(Category != "Total genes") %>%
  ggplot(aes(x = "", y = n_genes, fill = Category)) +
  geom_bar(stat = "identity", width=1, color="white") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "bottom") +
  # scale_fill_viridis(discrete = T) +
  geom_text(aes(label = paste0(Category, "\n", n_genes)), 
            position = position_stack(vjust = 0.5)) +
  labs(title = "Pangenome Summary Statistics",
       subtitle = "Number of genes in each category") +
  theme(legend.position = "bottom") +
  theme_cowplot(14)


# treemap

library(treemap)

pdf("exploration/PG_summary_treemap.pdf", width = 9, height = 7)
PG_summary %>% 
  filter(Category != "Total genes") %>%
  # filter(Category != "Cloud genes") %>% 
  treemap(index = c("Category", "n_genes"),
          vSize = "n_genes",
          type = "index",
          fontsize.labels = c(15, 12), # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
          fontcolor.labels = c("white", "black"), # Color of labels
          fontface.labels = c(2, 1), # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
          bg.labels = c("transparent"), # Background color of labels
          align.labels = list(
            c("center", "center"),
            c("right", "bottom")
          ), # Where to place labels in the rectangle?
          overlap.labels = 0.5, # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
          inflate.labels = F, # If true, labels are bigger when rectangle is bigger.
          palette = "Set1",
          title = "Mondriaan plot"
)

dev.off()


# ACC analysis  ---------------------------------------------------------------

# ACC means: 
# ACC data has been calculated in Python 
# with the script calculate_acc_parallel.py

acc_data = read_csv("data/panaroo_results/acc_data.csv")


acc_data.sum = acc_data %>% 
  group_by(n_genomes) %>% 
  summarise(mean_genes = mean(n_genes),
            sd_genes = sd(n_genes))

# Use log-log transformation to estimate initial values
log_fit <- lm(log(mean_genes) ~ log(n_genomes), data = acc_data.sum)

# Extract initial values from the linear model
initial_k <- exp(coef(log_fit)[1])
initial_b <- coef(log_fit)[2]

# Print initial values
initial_k
initial_b

# Fit Heap's law model using the improved initial values
heaps_model <- nls(mean_genes ~ k * n_genomes^b, data = acc_data.sum, 
                   start = list(k = initial_k, b = initial_b))

# Summarize the model to see the estimated parameters
summary(heaps_model)


# Extract the fitted parameters
k = coef(heaps_model)["k"]
b = coef(heaps_model)["b"]

# Create a sequence of genome numbers for plotting the fitted curve
genomes_seq = seq(min(acc_data.sum$n_genomes), max(acc_data.sum$n_genomes), length.out = 80)

# Calculate the fitted values using the estimated parameters
fitted_values = k * genomes_seq^b

# Add the fitted values to the original data frame
acc_data.sum = acc_data.sum %>%
  mutate(fitted_values = k * n_genomes^b)

# Plot the data
acc_data.sum %>%
  ggplot(aes(x = n_genomes, y = mean_genes)) +
  geom_line(aes(color = "Real Data")) +
  geom_line(data = data.frame(genomes_seq, fitted_values), 
            aes(x = genomes_seq, y = fitted_values, 
                color = "Fitted Data")) +
  geom_ribbon(aes(ymin = mean_genes - sd_genes, ymax = mean_genes + sd_genes), 
              fill = "grey", alpha = 0.5) +
  labs(title = "Heap's Law Fit to Pangenome Data",
       x = "Number of Genomes",
       y = "Number of Unique Genes",
       color = NULL) +
  theme_half_open() +
  background_grid() +
  scale_color_manual(values = c("Real Data" = "black", "Fitted Data" = "blue")) +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("Heap's law: y = ", round(k, 2), " * x^", 
                          round(b, 2)), 
           hjust = 1.1, vjust = 2, size = 5)

ggsave("exploration/Heap_law.pdf", width = 8, height = 7)



# Gene lists partition ----------------------------------------------------

gene_presence

core_genes = gene_presence %>% filter(gene_prop > 0.99) %>% pull(Gene)
shell_genes = gene_presence %>% filter(gene_prop <= 0.99, gene_prop > 0.15) %>% 
  pull(Gene)
cloud_genes = gene_presence %>% filter(gene_prop <= 0.15) %>% pull(Gene)

gene_presence = gene_presence %>% 
  mutate(gene_group = case_when(
    Gene %in% core_genes ~ "Core",
    Gene %in% shell_genes ~ "Shell",
    Gene %in% cloud_genes ~ "Cloud"
  )) 

gene_presence %>% 
  write_csv("tables/gene_groups.csv")


# gene partition by phylogroup --------------------------------------------

gene_pa_long = gene_pa %>% 
  pivot_longer(-Gene, names_to = "assembly_id_simp", values_to = "presence") %>% 
  filter(presence == 1) %>% 
  mutate(assembly_id_simp = case_when(str_detect(assembly_id_simp, "NT12") ~ str_remove(assembly_id_simp, "_.*"),
                                      TRUE ~ assembly_id_simp)) %>%
  left_join(metadata_final %>% select(assembly_id_simp, phylogroup)) 

gene_pa_long %>% 
  filter(is.na(phylogroup)) %>% 
  distinct(assembly_id_simp)


gene_count_phylo = gene_pa_long %>% 
  group_by(Gene, phylogroup) %>% 
  count()

phylo_counts = metadata_final %>% 
  count(phylogroup) %>% 
  rename(n_phylo = n)


gene_classes_phylogroup = gene_count_phylo %>% 
  arrange(desc(n)) %>% 
  left_join(phylo_counts) %>% 
  mutate(prop = n / n_phylo) %>% 
  mutate(class = case_when(
    prop > 0.95 ~ "Core",
    prop <= 0.95 & prop > 0.15 ~ "Shell",
    prop <= 0.15 ~ "Cloud"
  )) %>% 
  group_by(phylogroup, class) %>% 
  count() %>% 
  group_by(phylogroup) %>%
  mutate(prop = n / sum(n))

gene_classes_phylogroup %>% 
  write_csv("tables/gene_classes_phylogroup.csv")

# plot stacked barplot per phylogroup

gene_class_colors = c("Core" = "#158B81", 
                      "Cloud" = "#FF7F00", 
                      "Shell" = "#CD6090")

gene_classes_phylogroup %>% 
  ggplot(aes(phylogroup, n, fill = class)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Phylogroup", y = "Number of genes") +
  scale_fill_manual(values = gene_class_colors) +
  theme_cowplot(14) +
  geom_text(aes(label = scales::percent(prop)), 
            position = position_fill(vjust = 0.5))

ggsave("exploration/gene_classes_phylogroup.pdf", width = 8, height = 7)

gene_classes_phylogroup %>% 
  filter(class == "Core") %>%
  ggplot(aes(phylogroup, n, fill = phylogroup)) +
  geom_bar(stat = "identity") +
  labs(x = "Phylogroup", y = "Proportion of core genes") +
  scale_fill_manual(values = phylo_colors) 

gene_classes_phylogroup %>% 
  group_by(phylogroup) %>%
  mutate(total = sum(n)) %>%
  # filter(class == "Core") %>%
  ggplot(aes(x = n, y = total, color = phylogroup)) +
  geom_smooth(method = "lm", se = F, color = 'grey69') +
  geom_point(size = 4) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = phylo_colors) +
  facet_wrap(~class, scales = "free",
             ncol =1) +
  ggpubr::stat_cor(method = "pearson", size = 3,
                   label.y = 30000, label.x = 2500,
                   color = 'black') +
  labs(y = "Number of TOTAL gene families per phylogroup", 
       x = "Gene families in group") +
  theme_cowplot(14) 

ggsave("exploration/core_genes_phylogroup.pdf", width = 8, height = 7)

gene_classes_phylogroup %>% 
  group_by(phylogroup) %>%
  mutate(total = sum(n)) %>%
  # filter(class == "Core") %>%
  ggplot(aes(x = total, y = prop, color = phylogroup)) +
  geom_smooth(method = "lm", se = F, color = 'grey69') +
  geom_point(size = 4) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = phylo_colors) +
  facet_wrap(~class, scales = "free",
             ncol =1) +
  ggpubr::stat_cor(method = "pearson", size = 3,
                   label.y = 0.4, 
                   label.x = 28000,
                   color = 'black') +
  labs(x = "Total number of gene families per group", 
       y = "Proportion of gene families per group") +
  theme_cowplot(14) 

ggsave("exploration/core_genes_phylogroup_prop.pdf", width = 8, height = 7)



# read network output -----------------------------------------------------

## I'm reading here the paralogs to include them in the data analysis 

graph_data = read_csv("data/panaroo_results/final_graph.gml default node.csv")

graph_data = graph_data %>% 
  select(shared_name = `shared name`, everything(),
         -geneIDs, -genomeIDs)


graph_data %>% 
  count(paralog)

graph_data_seqs = graph_data %>% 
  # filter(!(shared_name %in% names(pangenome_ref))) %>% 
  select(shared_name, dna, protein)


graph_missing_seqs = graph_data_seqs %>%
  filter(!(shared_name %in% names(pangenome_ref))) 

graph_missing_seqs= graph_missing_seqs %>% 
  separate_rows(protein, sep = ";") %>% 
  separate_rows(dna, sep = ",") %>%
  group_by(shared_name) %>% 
  slice_head(n = 1) %>% 
  ungroup



# get protein sequences and save them as a fasta file
fasta = c()
for (i in 1:nrow(graph_missing_seqs)) {
  fasta = c(fasta, paste0(">", graph_missing_seqs$shared_name[i], "\n", 
                          graph_missing_seqs$protein[i]))
}

# save it as a fasta file with 
write(fasta, "data/panaroo_results/pan_genome_reference_aa_paralogs.fa")




# dna sequences

graph_missing_seqs = graph_missing_seqs %>% 
  separate_rows(dna, sep = ";") %>% 
  group_by(shared_name) %>%
  slice_head(n = 1) %>%
  ungroup() 

fasta = c()
for (i in 1:nrow(graph_missing_seqs)) {
  fasta = c(fasta, paste0(">", graph_missing_seqs$shared_name[i], "\n", 
                          graph_missing_seqs$dna[i]))
}

head(fasta)
# save it as a fasta file with
write(fasta, "data/panaroo_results/pan_genome_reference_paralogs.fa")


# remove vars to ƒree space

rm(fasta, graph_missing_seqs, graph_data_seqs, graph_data)
gc()


# GENE FUNCTIONS -----------------------------------------------------------

# note that the pangenome reference used from panaroo cluster paralogous genes 
# in only one entry, so the total number in the file is less than the total 
# number of gene families found in the gene presence absence matrix

# read fasta file
pangenome_ref = ape::read.FASTA("data/panaroo_results/pan_genome_reference_aa_TOTAL.fa")

str(pangenome_ref)

# number of genes in the pangenome reference
# 57219
refseqs_n = length(pangenome_ref)

### Proteinfer -------------------

# read the proteinfer results from main
proteinfer = read_tsv("data/proteinfer/proteinfer_output.tsv")

# read the proteinfer results from the paralogs
proteinfer_paralogs = read_tsv("data/proteinfer/proteinfer_paralogs_output.tsv")

proteinfer = proteinfer %>% 
  bind_rows(proteinfer_paralogs) 


proteinfer %>% 
  write_csv("tables/proteinfer_tidy_FINAL.csv")

proteinfer %>% 
  write_csv("GO_merge/proteinfer_tidy_FINAL.csv")


#### proteinfer metadata -----------

proteinfer_metadata = proteinfer %>% 
  distinct(predicted_label, description)

# how many genes with GO terms
proteinfer %>% 
  filter(str_detect(predicted_label, "GO")) %>% 
  distinct(sequence_name) %>% 
  count()

# 88.7% of genes had a GO term annotation
# 82208 / 92435 = 0.889

# how many genes with 1 GO term, with 2 GO terms, with 3... 
proteinfer %>% 
  filter(str_detect(predicted_label, "GO")) %>% 
  group_by(sequence_name) %>% 
  count() %>%
  ungroup %>% 
  count(n) %>% 
  rename(n_GO_terms = n,
         n_genes = nn) %>% 
  ggplot(aes(n_GO_terms, n_genes)) +
  geom_bar(stat = "identity") +
  labs(x = "Predicted GO terms", y = "Number of genes") +
  theme_cowplot(14)

ggsave("exploration/proteinfer/GO_terms_per_gene.pdf", width = 6, height = 4)


proteinfer %>% 
  filter(str_detect(predicted_label, "GO")) %>% 
  group_by(sequence_name) %>% 
  count() %>%
  ungroup %>% 
  count(n) %>% 
  rename(n_GO_terms = n,
         n_genes = nn) %>% 
  filter(n_GO_terms < 50) %>%
  ggplot(aes(n_GO_terms, n_genes)) +
  geom_bar(stat = "identity") +
  labs(x = "Predicted GO terms", y = "Number of genes") +
  theme_cowplot(14)

ggsave("exploration/proteinfer/GO_terms_per_gene_filtered.pdf", 
       width = 6, height = 4)
ggsave("exploration/proteinfer/GO_terms_per_gene_filtered.png", 
       width = 6, height = 4)


## GO terms distribution

meaningless_goterms = c("GO:0008150", "GO:0003674", "GO:0005575")

p1 = proteinfer %>% 
  filter(str_detect(predicted_label, "GO")) %>% 
  # filter(!predicted_label %in% meaningless_goterms) %>%
  count(predicted_label) %>% 
  arrange(desc(n)) %>% 
  head(20) %>% 
  ggplot(aes(reorder(predicted_label, n), n)) +
  geom_bar(stat = "identity", fill = c("#4166E0")) +
  coord_flip() +
  labs(x = "Predicted GO terms", y = "Number of genes") +
  theme_cowplot(14)

p2 = proteinfer %>%
  filter(str_detect(predicted_label, "GO")) %>% 
  # filter(!predicted_label %in% meaningless_goterms) %>%
  count(predicted_label, description) %>% 
  arrange(desc(n)) %>% 
  head(20) %>% 
  mutate(position = 1) %>% 
  ggplot(aes(x = position, y = fct_reorder(description, n))) +
  geom_text(aes(label = description), size = 3, hjust = 0) +
  xlim(0, 4) +
  # remove axis, labels and everything
  theme_void()

p1 + p2

ggsave("exploration/proteinfer/GO_terms_distribution_top20.pdf", 
       width = 8, height = 6)
ggsave("exploration/proteinfer/GO_terms_distribution_top20.png", 
       width = 8, height = 6)








### goPredSim --------------------

# read the goPredSim results
gopred = read_csv("data/goPredSim/goPredSim_final.csv") %>% 
  select(-`...1`)

# read the goPredSim results from the paralogs
gopred_paralogs = read_csv("data/goPredSim/goPredSim_paralogs.csv") %>% 
  select(-`...1`)


#### demultiplex from CD-HIT clusters -------------------------------------
# read the cd-hit clusters to demultiplex the data
cd_hit_clusters = read_table("data/goPredSim/cd_hit_clusters.tsv")
cd_hit_clusters_paralogs = read_table("data/goPredSim/cd_hit_paralogs_clusters.tsv")


# how many genes with GO terms
multicluster = cd_hit_clusters %>% 
  group_by(cluster) %>%
  count() %>% 
  filter(n > 1)

dim(multicluster)[1]

missing_gopred = cd_hit_clusters %>% 
  filter(cluster %in% multicluster$cluster) %>% 
  # filter(position == 0) %>% 
  rename(original_id = sequence) %>% 
  left_join(gopred) %>% 
  drop_na() %>% 
  select(-original_id, -position, -sequence_length) %>% 
  full_join(cd_hit_clusters %>% 
               filter(cluster %in% multicluster$cluster)) %>% 
  select(cluster, position, sequence, everything()) %>% 
  arrange(cluster) %>% 
  select(-cluster, -position) %>% 
  rename(original_id = sequence)

# join the missing values and the gopred values, remove duplicates
gopred = gopred %>% 
  bind_rows(missing_gopred) %>% 
  distinct(original_id, .keep_all = TRUE)


# demultiplex paralogs now
multicluster_paralogs = cd_hit_clusters_paralogs %>% 
  group_by(cluster) %>%
  count() %>% 
  filter(n > 1)

missing_gopred_paralogs = cd_hit_clusters_paralogs %>%
  filter(cluster %in% multicluster_paralogs$cluster) %>% 
  # filter(position == 0) %>% 
  rename(original_id = sequence) %>% 
  left_join(gopred_paralogs) %>% 
  drop_na() %>% 
  select(-original_id, -position, -sequence_length) %>% 
  full_join(cd_hit_clusters_paralogs %>% 
               filter(cluster %in% multicluster_paralogs$cluster)) %>% 
  select(cluster, position, sequence, everything()) %>% 
  arrange(cluster) %>% 
  select(-cluster, -position) %>% 
  rename(original_id = sequence)


gopred_paralogs = gopred_paralogs %>%
  bind_rows(missing_gopred_paralogs) %>% 
  distinct(original_id, .keep_all = TRUE)


## JOIN THE TWO DATASETS
gopred = gopred %>% 
  bind_rows(gopred_paralogs) %>% 
  distinct(original_id, .keep_all = TRUE)


gopred %>% 
  write_csv("tables/gopredsim_tidy_FINAL.csv")
gopred %>% 
  write_csv("GO_merge/gopredsim_tidy_FINAL.csv")


# distance density plot
gopred %>% 
  select(original_id, k_nn_1_distance, k_nn_2_distance, k_nn_3_distance) %>% 
  pivot_longer(-original_id, names_to = "k_nn", values_to = "distance") %>%
  ggplot(aes(distance, color = k_nn, fill = k_nn)) +
  geom_density(alpha = .3) +
  labs(x = "Euclidean distance", y = "Density") +
  theme_cowplot(14)

ggsave("exploration/gopredsim/gopredsim_distances.pdf", width = 6, height = 4)



### 


gopred_dist_sum = gopred %>% 
  select(original_id, k_nn_1_distance, k_nn_2_distance, k_nn_3_distance) %>% 
  pivot_longer(-original_id, names_to = "k_nn", values_to = "distance") %>%
  # filter(distance < 2) %>% 
  group_by(k_nn) %>% 
  summarise(Mean = mean(distance),
            SD = sd(distance),
            Median = median(distance),
            Min = min(distance),
            Max = max(distance))


gopred %>% 
  select(original_id, k_nn_1_distance, k_nn_2_distance, k_nn_3_distance) %>% 
  pivot_longer(-original_id, names_to = "k_nn", values_to = "distance") %>%
  filter(distance < 2) %>% 
  ggplot(aes(distance, color = k_nn, fill = k_nn)) +
  geom_density(alpha = .3) +
  # vertical lines as median from gopred_dist_sum
  geom_vline(data = gopred_dist_sum, 
             aes(xintercept = Median, color = k_nn), linetype = "dashed") +
  geom_vline(xintercept = 1.05) +
  labs(x = "Euclidean distance", y = "Density") +
  theme_cowplot(14)

ggsave("exploration/gopredsim/gopredsim_distances_capped.pdf", width = 6, height = 4)
ggsave("exploration/gopredsim/gopredsim_distances_capped.png", width = 6, height = 4)
# make the list longer

gopred_dist_long = gopred %>% 
  select(original_id, k_nn_1_distance, k_nn_2_distance, k_nn_3_distance) %>%
  pivot_longer(-original_id, names_to = "k_nn", values_to = "distance") %>% 
  mutate(kNN = case_when(
    k_nn == "k_nn_1_distance" ~ 1,
    k_nn == "k_nn_2_distance" ~ 2,
    k_nn == "k_nn_3_distance" ~ 3
  ))


gopred_annotations_long = gopred %>% 
  select(original_id, k_nn_1_annotations, k_nn_2_annotations, k_nn_3_annotations) %>%
  pivot_longer(-original_id, names_to = "k_nn", values_to = "annotations") %>% 
  mutate(kNN = case_when(
    k_nn == "k_nn_1_annotations" ~ 1,
    k_nn == "k_nn_2_annotations" ~ 2,
    k_nn == "k_nn_3_annotations" ~ 3
  ))

gopred_long = gopred_dist_long %>%
  select(-k_nn) %>% 
  left_join(gopred_annotations_long) %>% 
  separate_rows(annotations, sep = ";")



gopred_long %>%
  filter(distance < 1.05) %>% 
  distinct(original_id) %>% 
  count()



### Interproscan -------------------------

# read the interproscan results
interpro = read_tsv("data/interpro/pan_genome_reference_aa.fa.tsv",
                    col_names = FALSE) %>% 
  rename(original_id = X1,
         protein_length = X3, 
         database = X4,
         database_id = X5,
         domain = X6,
         start_location = X7,
         stop_location = X8,
         evalue = X9, 
         go_term = X14,
         pathway = X15) %>%
  select(
    original_id, protein_length, database, database_id,
    evalue, domain, go_term,
    start_location, stop_location, pathway
  )

# read the interproscan results from the paralogs
interpro_paralogs = read_tsv("data/interpro/pan_genome_reference_aa_paralogs.fa.tsv",
                    col_names = FALSE) %>% 
  rename(original_id = X1,
         protein_length = X3, 
         database = X4,
         database_id = X5,
         domain = X6,
         start_location = X7,
         stop_location = X8,
         evalue = X9, 
         go_term = X14,
         pathway = X15) %>%
  select(
    original_id, protein_length, database, database_id,
    evalue, domain, go_term,
    start_location, stop_location, pathway
  )


# merge both datasets
interpro = interpro %>% 
  bind_rows(interpro_paralogs) 


interpro %>% 
  write_csv("tables/interproscan_RAW.csv")


# tidy the interpro data


interpro_tidy = interpro %>% 
  filter(go_term != "-") %>%
  mutate(evalue = as.numeric(evalue)) %>%
  separate_rows(go_term, sep = "\\|") %>% 
  mutate(go_term = str_remove(go_term, "\\(.*\\)")) %>% 
  group_by(original_id, protein_length, go_term) %>%
  drop_na(evalue) %>% 
  slice_min(evalue) %>%
  ungroup %>% 
  filter(evalue < 1e-5)


interpro_tidy %>%
  distinct(original_id) %>%
  count()


interpro_tidy %>%
  write_csv("tables/interproscan_tidy_FINAL.csv")
interpro_tidy %>%
  write_csv("GO_merge/interproscan_tidy_FINAL.csv")









### eggnog -----------------------


# read the eggnog results
eggnog = read_delim("data/eggnog/complete_PG_eggnogmap.emapper.annotations", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE, skip = 4) %>% 
  rename(original_id = `#query`)


# read eggnog from paralogs
eggnog_paralogs = read_delim("data/eggnog/complete_PG_eggnogmap_paralogs.emapper.annotations", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE, skip = 4) %>% 
  rename(original_id = `#query`)


# merge both datasets
eggnog = eggnog %>% 
  bind_rows(eggnog_paralogs) 

# this dataset is IMPORTANT
# remove unmapped COG terms from the dataset
eggnog_COG = eggnog %>% 
  filter(COG_category != "-")

# remove unmapped GO terms from dataset
eggnog = eggnog %>% 
  filter(GOs != "-")  %>% 
  filter(evalue < 1e-5)


eggnog %>% 
  distinct(original_id) %>% 
  count()


eggnog %>% 
  write_csv("tables/eggnog_tidy_FINAL.csv")
eggnog %>% 
  write_csv("GO_merge/eggnog_tidy_FINAL.csv")


# # # # # # # # 

eggnog %>% 
  filter(GOs != "-") %>% 
  distinct(original_id, .keep_all = T) %>% 
  count()


eggnog %>% 
  filter(GOs == "-")

eggnog %>% 
  filter(Preferred_name == "-") %>% 
  filter(!str_detect(original_id, "group"))



### joint analysis -------------------

#### number of gene families with GO terms -------------

prot_n = proteinfer %>% 
  filter(str_detect(predicted_label, "GO")) %>% 
  filter(confidence > 0.4) %>% 
  distinct(sequence_name) %>% 
  count() %>% pull(n)

gopred_n = gopred_long %>%
  filter(distance < 1.05) %>% 
  distinct(original_id) %>% 
  count() %>% pull(n)


interpro_n = interpro_tidy %>%
  distinct(original_id) %>% 
  count() %>% pull(n)


eggnog_n = eggnog %>%
  filter(GOs != "-") %>% 
  distinct(original_id) %>% 
  count() %>% pull(n)


# create a tibble with results
total_n = length(pangenome_ref)
gene_families_GO = tibble(
  method = c("Proteinfer", "goPredSim", "Interproscan", "Eggnog"),
  n_genes = c(prot_n, gopred_n, interpro_n, eggnog_n)
)

# set 1
colors_methods = c("Proteinfer" = "#1f77b4", 
            "eggNOG" =  "#ff7f0e", 
            "goPredSim" = "#2ca02c",
            "Interproscan" = "#d62728") 

gene_families_GO %>% 
  mutate(n_prop = n_genes / total_n) %>% 
  ggplot(aes(fct_reorder(method, n_prop, .desc = TRUE), n_prop, 
             fill = method)) +
  geom_bar(stat = "identity", show.legend = F) +
  # annotate with the prop number on each barplot
  geom_text(aes(label = scales::percent(n_prop)), 
            position = position_stack(vjust = 0.5)) +

  labs(x = "Method", y = "Proportion of genes with GO terms") +
  scale_fill_manual(values = colors_methods ) +
  theme_half_open(14) +
  ylim(0, 1)

ggsave("exploration/gene_families_GO_methods.pdf", width = 6, height = 4)
ggsave("exploration/gene_families_GO_methods.png", width = 6, height = 4)

#### venn diagram of genes in each method -------------------

# create a list of genes with GO terms
prot_genes = proteinfer %>% 
  filter(str_detect(predicted_label, "GO")) %>% 
  distinct(sequence_name) %>% 
  pull(sequence_name)


gopred_genes = gopred_long %>%
  filter(distance < 1.05) %>% 
  distinct(original_id) %>% 
  pull(original_id)


interpro_genes = interpro_tidy %>%
  distinct(original_id) %>% 
  pull(original_id)


eggnog_genes = eggnog %>%
  filter(GOs != "-") %>% 
  distinct(original_id) %>% 
  pull(original_id)


# create the venn diagram
# library(UpSetR)
library(ComplexUpset)

# create a tibble with genes, methods as columns, and 1 and 0 as values
genes_methods = tibble(
  original_id = names(pangenome_ref) ,
  Proteinfer = as.integer(names(pangenome_ref) %in% prot_genes),
  goPredSim = as.integer(names(pangenome_ref) %in% gopred_genes),
  Interproscan = as.integer(names(pangenome_ref) %in% interpro_genes),
  Eggnog = as.integer(names(pangenome_ref) %in% eggnog_genes)
)

# make a list of the results
genes_methods_list = list(
  Proteinfer = prot_genes,
  goPredSim = gopred_genes,
  Interproscan = interpro_genes,
  Eggnog = eggnog_genes
)
 
upset(fromList(genes_methods_list),
        order.by = "freq", 
        decreasing = T, 
        mb.ratio = c(0.6, 0.4),
        number.angles = 0, 
        text.scale = 1.1, 
        point.size = 2.8, 
        line.size = 1
  )


genes_methods_classes = genes_methods %>% 
  left_join(gene_presence %>% 
              rename(original_id = Gene)) 

method_names = c("Proteinfer", "goPredSim", "Interproscan", "Eggnog")

upset(
  genes_methods_classes,
  method_names,
  base_annotations = list(
    'Intersection size'= intersection_size(
      counts = FALSE,
      mapping = aes(fill = gene_group)
    ) 
  ),
  width_ratio=0.1
)



gene_fun_upset = function(mode = "exclusive_intersection") upset(
  genes_methods_classes, 
  method_names, 
  mode = mode, 
  set_sizes = FALSE,
  encode_sets = FALSE,
  base_annotations=list(
    'Size'=(
      intersection_size(
        mode = mode,
        mapping = aes(fill=exclusive_intersection),
        size = 0,
        text = list(check_overlap=TRUE)
      ) + scale_fill_venn_mix(
        data = genes_methods_classes,
        guide = 'none',
        colors = gene_class_colors
      )
    )
  )
)

gene_fun_upset("exclusive_intersection")


upset(
  genes_methods_classes,
  method_names,
  # mode = "inclusive_union",
  base_annotations = list(
    'Pangenome class' = intersection_size(
      counts = T,
      mapping = aes(fill = gene_group),
        text=list(
          vjust=-0.1,
          hjust=-0.1,
          angle = 45,
          size = 3
        ) ) +
      scale_fill_manual(values = gene_class_colors)
    # ,
    # 'Intersection ratio' = intersection_ratio(
    #   mapping = aes(fill = gene_group),
    #   text=list(
    #     vjust=-0.,
    #     hjust=-0.,
    #     angle = 90,
    #     size = 3
    #   ) 
    # ) +
    #   scale_fill_manual(values = gene_class_colors)
  ),
  width_ratio=0.1
)

# save plot
ggsave("exploration/gene_families_GO_venn.pdf", width = 12, height = 8)
ggsave("exploration/gene_families_GO_venn.png", width = 15, height = 10)




upset(
  genes_methods_classes,
  method_names,
  # mode = "inclusive_union",
  base_annotations = list(
    'Pangenome class' = intersection_size(
      counts = T,
      mapping = aes(fill = gene_group),
      text=list(
        vjust=-0.1,
        hjust=-0.1,
        angle = 45,
        size = 2
      ),
      text_mapping=aes(label=paste0(
        !!upset_text_percentage(),
        '\n',
        '(',
        !!get_size_mode('exclusive_intersection'),
        ')'
      ))) +
      scale_fill_manual(values = gene_class_colors)
  ),
  width_ratio=0.1
)




x# save the gene functions data
genes_methods_classes %>% 
  write_csv("tables/function_methods_data.csv")



stats_upset = upset_test(genes_methods_classes,
           method_names)




# IC content analysis -----------------------------------------------------

# the data from the IC content comes from the folder GO_merge. There
# I have calcualted the IC content for each gene GO term

ic_content = read_csv("GO_merge/IC_GOterms.csv") %>% 
  rename(GO = Query_GO)

all_tidy_func = read_csv("GO_merge/all_tidy_FINAL.csv") %>% 
  select(gene_family, GO, method)

gene_presence


ic_tbl = all_tidy_func %>% 
  left_join(ic_content) %>% 
  left_join(gene_presence %>% 
              rename(gene_family = Gene)) %>% 
  drop_na(IC, gene_group)


# some summary statistics

ic_sum_stats = ic_tbl %>%
  group_by(method, gene_group) %>%
  summarise(
    mean_IC = mean(IC),
    sd_IC = sd(IC),
    SEM_IC = sd(IC) / sqrt(n()),
    median_IC = median(IC),
    MAD_IC = mad(IC),
    min_IC = min(IC),
    max_IC = max(IC),
    # Calculate the confidence interval for the median
    median_CI = list(MedianCI(IC, conf.level = 0.95)), # Choose confidence level
    # Calculate the standard error of the median 
    SEMd_IC = abs(MedianCI(IC, conf.level = 0.95)[2] - 
                    MedianCI(IC, conf.level = 0.95)[1]) / (2 * qnorm(0.975))
  ) %>%
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud"))) %>% 
  ungroup





# check if data is normally distributed
sample(ic_tbl$IC, 1000) %>% shapiro.test(.)
# the data IS NOT normally distributed

ic_wilcox_general_stats = ic_tbl %>%
  group_by(method) %>% 
  mutate(gene_group = factor(gene_group, levels = c("Core", "Shell", "Cloud"))) %>%
  # Mann-Whitney U test 
  pairwise_wilcox_test(IC ~ gene_group, p.adjust.method = "fdr") %>% 
  add_xy_position()

ic_wilcox_general_stats_2 = ic_tbl %>%
  group_by(gene_group) %>% 
  mutate(gene_group = factor(gene_group, levels = c("Core", "Shell", "Cloud"))) %>%
  # Mann-Whitney U test 
  pairwise_wilcox_test(IC ~ method, p.adjust.method = "fdr") 


list_of_datasets = list(
  "wilcox per method" = ic_wilcox_general_stats,
  "wilcox per group" = ic_wilcox_general_stats_2

)


write.xlsx(list_of_datasets, "tables/IC_content_general_stats.xlsx")



# with mean and SEM  
ic_sum_stats %>% 
  ggplot(aes(method, mean_IC, fill = gene_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_IC - SEM_IC, ymax = mean_IC + SEM_IC),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(x = "Method", y = "Overall mean IC content ± SEM") +
  theme_cowplot(14) +
  scale_fill_manual(values = gene_class_colors)

ggsave("exploration/gene_functions_analysis/IC_content_mean.pdf", 
       width = 8, height = 6)


# median and SEMd_IC
ic_sum_stats %>% 
  ggplot(aes(method, median_IC, fill = gene_group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = median_IC - SEMd_IC, ymax = median_IC + SEMd_IC),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(x = "Method", y = "Overall median IC content ± SEMd") +
  theme_cowplot(14) +
  scale_fill_manual(values = gene_class_colors)

ggsave("exploration/gene_functions_analysis/IC_content_median.pdf", 
       width = 8, height = 6)


### boxplots -----------

ic_tbl %>% 
  ggplot(aes(gene_group, IC, fill = gene_group)) +
  geom_boxplot() +
  labs(x = "Method", y = "Information content") +
  theme_cowplot(14) +
  scale_fill_manual(values = gene_class_colors) +
  facet_wrap(~method, scales = "free_x", ncol = 4) +
  # remove x-axis names
  background_grid() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  panel_border()

ggsave("exploration/gene_functions_analysis/IC_content_boxplot.pdf", 
       width = 8, height = 6)
ggsave("exploration/gene_functions_analysis/IC_content_boxplot.png", 
       width = 8, height = 6)


ggpubr::ggboxplot(ic_tbl, x = "gene_group", y = "IC",
                  color = "black", 
                  fill = "gene_group",
                  facet.by = "method",
                  # facet in 4 columns
                  ncol = 4,
                  palette = gene_class_colors) +
  ggpubr::stat_pvalue_manual(ic_wilcox_general_stats, 
                             label = "p.adj.signif",
                             tip.length = 0.01)


ggsave("exploration/gene_functions_analysis/IC_content_boxplot_facet.pdf",
       height = 8, width = 12)
ggsave("exploration/gene_functions_analysis/IC_content_boxplot_facet.png",
       height = 8, width = 12)






## Max IC per gene analysis --------------------------------------------

# get the max IC per gene
max_ic = ic_tbl %>% 
  group_by(gene_family, method) %>% 
  summarise(max_IC = max(IC)) %>% 
  ungroup() %>% 
  left_join(gene_presence %>% 
              rename(gene_family = Gene)) %>% 
  drop_na(max_IC, gene_group)


# some summary statistics

max_ic_sum_stats = max_ic %>%
  group_by(method, gene_group) %>%
  summarise(
    mean_IC = mean(max_IC),
    sd_IC = sd(max_IC),
    SEM_IC = sd(max_IC) / sqrt(n()),
    median_IC = median(max_IC),
    MAD_IC = mad(max_IC),
    min_IC = min(max_IC),
    max_IC = max(max_IC)) %>%
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud"))) %>% 
  ungroup

# Wilcoxon test
max_ic_wilcox_general_stats = max_ic %>%
  group_by(method) %>% 
  mutate(gene_group = factor(gene_group, levels = c("Core", "Shell", "Cloud"))) %>%
  # Mann-Whitney U test 
  pairwise_wilcox_test(max_IC ~ gene_group, p.adjust.method = "fdr") %>% 
  add_xy_position()

max_ic_wilcox_general_stats_2 = max_ic %>%
  group_by(gene_group) %>% 
  mutate(gene_group = factor(gene_group, levels = c("Core", "Shell", "Cloud"))) %>%
  # Mann-Whitney U test 
  pairwise_wilcox_test(max_IC ~ method, p.adjust.method = "fdr")

list_of_datasets = list(
  "wilcox per method" = max_ic_wilcox_general_stats,
  "wilcox per group" = max_ic_wilcox_general_stats_2
)

write.xlsx(list_of_datasets, "tables/IC_content_maxIC_stats.xlsx")



# boxplots

max_ic %>% 
  ggplot(aes(gene_group, max_IC, fill = gene_group)) +
  geom_boxplot() +
  labs(x = "Method", y = "Max IC content") +
  theme_cowplot(14) +
  scale_fill_manual(values = gene_class_colors) +
  facet_wrap(~method, scales = "free_x", ncol = 4) 


ggpubr::ggboxplot(max_ic, x = "gene_group", y = "max_IC",
                  color = "black", 
                  fill = "gene_group",
                  facet.by = "method",
                  # facet in 4 columns
                  ncol = 4,
                  palette = gene_class_colors) +
  ggpubr::stat_pvalue_manual(max_ic_wilcox_general_stats, 
                             label = "p.adj.signif",
                             tip.length = 0.01)


ggsave("exploration/gene_functions_analysis/IC_content_maxIC_boxplot.pdf",
       height = 8, width = 12)
ggsave("exploration/gene_functions_analysis/IC_content_maxIC_boxplot.png",
       height = 8, width = 12)


## analysis of shared genes -----------------------------------

genes_shared = max_ic %>% 
  group_by(gene_family) %>% 
  count() %>% 
  filter(n == 4) %>%
  pull(gene_family) 


max_ic_shared = max_ic %>% 
  filter(gene_family %in% genes_shared) 

















# T5 UMAP projections -----------------------------------------------------


umap_proj = read_csv("data/t5_embeddings/umap_projections/projected_embeddings_file.csv") %>% 
  rename(id = `...1`)


umap_proj %>% 
  mutate(log_seq_len = log10(sequence_length)) %>% 
  ggplot(aes(x = component_0, component_1,
             color = log_seq_len)) +
  geom_point(alpha = 0.2) +
  scale_color_viridis_c("Log10 \n(sequence length)") +
  theme_cowplot(14) +
  labs(title = "UMAP projection of T5 embeddings",
       subtitle = "Color by log10 of sequence length",
       x = "Component 0", y = "Component 1")

ggsave("exploration/t5_embeddings/umap_projection_T5.pdf", 
       width = 8, height = 6)
ggsave("exploration/t5_embeddings/umap_projection_T5.png", 
       width = 8, height = 6)


umap_proj_cloud = umap_proj %>% 
  rename(Gene = original_id) %>% 
  left_join(gene_presence) %>% 
  filter(gene_group == "Cloud")

umap_proj_core = umap_proj %>% 
  rename(Gene = original_id) %>% 
  left_join(gene_presence) %>% 
  filter(gene_group == "Core")

umap_proj_shell = umap_proj %>%
  rename(Gene = original_id) %>% 
  left_join(gene_presence) %>% 
  filter(gene_group == "Shell")


umap_proj_cloud %>% 
  ggplot(aes(x = component_0, component_1,
             color = gene_group)) +
  geom_point(alpha = 0.2) +
  geom_point(data = umap_proj_shell, alpha = 0.4) +
  geom_point(data = umap_proj_core, alpha = 0.4) +
  scale_color_manual(values = gene_class_colors) +
  theme_cowplot(14) +
  labs(title = "UMAP projection of T5 embeddings",
       subtitle = "Color by gene group",
       x = "Component 0", y = "Component 1")

ggsave("exploration/t5_embeddings/umap_projection_T5_gene_groups.pdf", 
       width = 8, height = 6)
ggsave("exploration/t5_embeddings/umap_projection_T5_gene_groups.png", 
       width = 8, height = 6)


umap_proj_cloud %>% 
  ggplot(aes(x = component_0, component_1,
             color = gene_prop)) +
  geom_point(alpha = 0.2) +
  geom_point(data = umap_proj_shell, alpha = 0.4) +
  geom_point(data = umap_proj_core, alpha = 0.4) +
  scale_color_viridis_c("Gene proportion") +
  theme_cowplot(14)

ggsave("exploration/t5_embeddings/umap_projection_T5_gene_prop.pdf", 
       width = 8, height = 6)
ggsave("exploration/t5_embeddings/umap_projection_T5_gene_prop.png",
       width = 8, height = 6)


# select the gene families that appear only once

unique_genes = gene_pa_long %>% 
  count(Gene) %>% 
  filter(n == 1)


unique_genes = gene_pa_long %>% 
  filter(Gene %in% unique_genes$Gene) %>% 
  select(Gene, phylogroup)


unique_genes %>%
  rename(original_id = Gene) %>% 
  left_join(umap_proj) %>% 
  drop_na(component_0) %>% 
  ggplot(aes(x = component_0, component_1,
             color = phylogroup)) +
  geom_point(alpha = 0.7,
             size = 1) +
  scale_color_manual(values = phylo_colors) +
  theme_cowplot(14) +
  labs(title = "UMAP projection of T5 embeddings",
       subtitle = "Unique gene families",
       x = "Component 0", y = "Component 1")

ggsave("exploration/t5_embeddings/umap_projection_T5_unique_genes.pdf", 
       width = 8, height = 6)
ggsave("exploration/t5_embeddings/umap_projection_T5_unique_genes.png",
       width = 8, height = 6)


unique_genes %>%
  rename(original_id = Gene) %>% 
  left_join(umap_proj) %>% 
  drop_na(component_0) %>% 
  ggplot(aes(x = component_0, component_1,
             color = phylogroup)) +
  geom_point(alpha = 0.7,
             size = 1) +
  scale_color_manual(values = phylo_colors) +
  theme_cowplot(14) +
  labs(title = "UMAP projection of T5 embeddings",
       subtitle = "Unique gene families",
       x = "Component 0", y = "Component 1") +
  facet_wrap(~phylogroup)

ggsave("exploration/t5_embeddings/umap_projection_T5_unique_genes_facet.pdf",
       width = 12, height = 10)
ggsave("exploration/t5_embeddings/umap_projection_T5_unique_genes_facet.png",
       width = 12, height = 10)



# fix T5 embeddings file --------------------------------------------------

# read the embeddings
t5_embeddings = 
  read_csv("data/t5_embeddings/t5_embeddings/reduced_embeddings_file.csv") %>% 
  rename(id = `...1`)

t5_embeddings = t5_embeddings %>% 
  left_join(umap_proj %>% select(id, sequence_length, original_id)) %>% 
  select(original_id, everything())

missing_t5_emb = pca_embeddings_clusters %>% 
  filter(cluster %in% multicluster$cluster) %>% 
  left_join(t5_embeddings) %>% 
  drop_na() %>% 
  select(-original_id, -position, -sequence_length) %>% 
  full_join(pca_embeddings_clusters %>% 
              filter(cluster %in% multicluster$cluster)) %>% 
  select(cluster, position, original_id, everything()) %>% 
  arrange(cluster) %>% 
  drop_na(id) %>% 
  select(-cluster, -position, -id) %>% 
  distinct(original_id, .keep_all = TRUE)

t5_embeddings = t5_embeddings %>%
  bind_rows(missing_t5_emb) %>% 
  distinct(original_id, .keep_all = TRUE)

t5_embeddings = t5_embeddings %>% 
  left_join(gene_presence %>% 
              rename(original_id =Gene)) 

t5_embeddings = t5_embeddings %>% 
  select(original_id, gene_group, everything()) %>% 
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud"))) %>%
  arrange(desc(gene_group)) %>% 
  select(-id)


# save embeddings in csv
t5_embeddings %>% 
  select(-sequence_length, -gene_presence, -gene_prop) %>% 
  write_csv("data/t5_embeddings/t5_embeddings/t5_embeddings_METADATA_FINAL.csv")


# PCA projections -------------------------------------------------------

pca_embeddings = read_csv("exploration/t5_embeddings/PCA/pca_embeddings.csv") %>% 
  rename(original_id = Gene,
         pc1 = x, 
         pc2 = y) %>% 
  drop_na(gene_group)

per_var = c("9.28%", "4.84%")

### OLD DATA VERSION, DO NOT USE
### THE NEW ONE HAS BEEN SCALED BEFORE PERFOMING THE PCA
# pca_embeddings = read_csv("data/t5_embeddings/t5_embeddings/pca_embeddings.csv") %>% 
#   rename(
#          pc1 = `principal component 1`,
#          pc2 = `principal component 2`)
# 
# pca_embeddings = pca_embeddings %>% 
#   left_join(umap_proj %>% select(id, sequence_length, original_id)) %>% 
#   select(pc1, pc2, original_id)
# 
# 
# #### IMPORTANT
# ####
# ####
# 
# # percentage of variation explained by the two PCs: [0.12782702, 0.05466168]
# 
# # load clusters and demultiplex
# pca_embeddings_clusters = read_table("data/t5_embeddings/t5_embeddings/PCA_embeddings_clusters.tsv")
# 
# multicluster = pca_embeddings_clusters %>% 
#   group_by(cluster) %>%
#   count() %>%
#   filter(n > 1)
# 
# dim(multicluster)[1]
# 
# missing_pca_emb = pca_embeddings_clusters %>% 
#   filter(cluster %in% multicluster$cluster) %>% 
#   left_join(pca_embeddings) %>% 
#   drop_na() %>% 
#   select(-original_id, -position, -sequence_length) %>% 
#   full_join(pca_embeddings_clusters %>% 
#               filter(cluster %in% multicluster$cluster)) %>% 
#   select(cluster, position, original_id, everything()) %>% 
#   arrange(cluster) %>% 
#   drop_na(id) %>% 
#   select(-cluster, -position, -`...1`, -id) %>% 
#   distinct(original_id, .keep_all = TRUE)
# 
# missing_pca_emb
# 
# pca_embeddings = pca_embeddings %>%
#   bind_rows(missing_pca_emb) %>% 
#   distinct(original_id, .keep_all = TRUE)
# 
# pca_embeddings = pca_embeddings %>% 
#   left_join(gene_presence %>% 
#               rename(original_id =Gene)) 
# 
# pca_embeddings = pca_embeddings %>% 
#   mutate(gene_group = factor(gene_group, 
#                              levels = c("Core", "Shell", "Cloud"))) %>%
#   arrange(desc(gene_group))
# 
# pca_embeddings = pca_embeddings %>% 
#   drop_na(gene_group)
# 
# # write the pca embeddings to a file
# 
# pca_embeddings %>% 
#   write_csv("tables/pca_embeddings.csv")

#### general PCA projection --------------

pca_embeddings %>% 
  ggplot(aes(x = pc1, pc2,
             color = gene_group)) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = gene_class_colors) +
  theme_cowplot(14) +
  labs(title = "PCA projection of T5 embeddings",
       subtitle = "Color by gene group",
       x = glue("Principal component 1 ({per_var[1]})"), 
       y = glue("Principal component 2 ({per_var[2]})"))

ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups.pdf",
       width = 8, height = 6)
ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups.png",
       width = 8, height = 6)

#### PCA projection + 2D density ---------

p = pca_embeddings %>% 
  ggplot(aes(x = pc1, y = pc2,
             color = gene_group)) +
  geom_point(alpha = 0.2) +
  geom_density_2d(contour_var = "ndensity",
                  linewidth = 0.35,
                  alpha = 0.4) +
  scale_color_manual(values = gene_class_colors) +
  theme_cowplot(14) +
  labs(title = "PCA projection of T5 embeddings",
       subtitle = "Color by gene group",
       x = glue("Principal component 1 ({per_var[1]})"), 
       y = glue("Principal component 2 ({per_var[2]})"))

p = ggMarginal(p,
           type = "density",
           groupColour = TRUE,
           groupFill = TRUE)

p

ggsave(plot = p, "exploration/t5_embeddings/pca_projection_T5_gene_groups.pdf",
       width = 8, height = 6)
ggsave(plot = p, "exploration/t5_embeddings/pca_projection_T5_gene_groups.png",
       width = 8, height = 6)


#### PCA facet -----------------------------------

pca_embeddings %>% 
  ggplot(aes(x = pc1, pc2,
             color = gene_group)) +
  geom_point(alpha = 0.1) +
  scale_color_manual(values = gene_class_colors) +
  theme_cowplot(14) +
  labs(title = "PCA projection of T5 embeddings",
       subtitle = "Color by gene group",
       x = glue("Principal component 1 ({per_var[1]})"), 
       y = glue("Principal component 2 ({per_var[2]})")) +
  facet_wrap(~gene_group)

ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_facet.pdf",
       width = 8, height = 6)
ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_facet.png",
       width = 8, height = 6)

# 2D density plot
pca_embeddings %>% 
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud"))) %>% 
  ggplot(aes(x = pc1, y = pc2,
             color = gene_group)) +
  geom_density_2d(contour_var = "ndensity",
                  linewidth = 0.35) +
  theme_cowplot(14) +
  scale_color_manual(values = gene_class_colors) +
  labs(title = "PCA projection of T5 embeddings",
       subtitle = "Density plot of Cloud genes",
       x = glue("Principal component 1 ({per_var[1]})"), 
       y = glue("Principal component 2 ({per_var[2]})")) +
  facet_wrap(~gene_group)

ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_density.pdf",
       width = 8, height = 6)
ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_density.png",
       width = 8, height = 6)

# use geom_bin2d

pca_embeddings %>% 
  ggplot(aes(x = pc1, y = pc2)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", 
                  contour = T) +
  theme_cowplot(14) +
  labs(title = "PCA projection of T5 embeddings",
       subtitle = "Density plot of Cloud genes",
       x = glue("Principal component 1 ({per_var[1]})"), 
       y = glue("Principal component 2 ({per_var[2]})")) +
  facet_wrap(~gene_group)

ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_density2.pdf",
       width = 8, height = 6)
ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_density2.png",
       width = 8, height = 6)


### PCA axes boxplots ----------------------------------------------

# calculate stats for the PCA axes
pca_stats = pca_embeddings %>% 
  pivot_longer(-c(original_id, gene_presence, gene_prop, gene_group),
               names_to = "PC",
               values_to = "value") %>%
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud")),
         PC = factor(PC)) %>% 
  group_by(PC) %>% 
  rstatix::wilcox_test(value ~ gene_group) %>% 
  add_xy_position(x = "PC")


pca_embeddings %>% 
  pivot_longer(-c(original_id, gene_presence, gene_prop, gene_group),
               names_to = "PC",
               values_to = "value") %>%
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud")),
         PC = factor(PC)) %>% 
  filter(value > quantile(value, 0.025)) %>%
  ggplot(aes(x = PC, y = value,
             fill = gene_group)) +
  # geom_boxplot() 
  geom_violin() +
  geom_boxplot(width = 0.07, 
               position = position_dodge(0.9)) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)), color = "none") +
  scale_fill_manual(values = gene_class_colors) +
  labs(
    y = "Principal Component value",
    x = NULL,
    fill = "Group",
    caption = "95% CI"
  )

# ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_axes.pdf",
#        width = 8, height = 6)
# ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_axes.png",
#        width = 8, height = 6)


# ggpubr ggviolin

pca_embeddings %>% 
  pivot_longer(-c(original_id, gene_presence, gene_prop, gene_group),
               names_to = "PC",
               values_to = "value") %>%
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud")),
         PC = factor(PC)) %>% 
  filter(value > quantile(value, 0.025)) %>%
  ggviolin(x = "PC", y = "value",
            color = "black", 
            fill = "gene_group",
           add = "boxplot",
           add.params = list(
                             width = 0.1),
            # facet.by = "PC",
            # facet in 4 columns
            ncol = 2,
            palette = gene_class_colors) +
  ggpubr::stat_pvalue_manual(pca_stats, 
                             label = "p.adj.signif",
                             tip.length = 0.01) +
  labs(
    y = "Principal Component value",
    x = NULL,
    fill = "Group"
  )

ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_axes.pdf",
       width = 7, height = 5)
ggsave("exploration/t5_embeddings/pca_projection_T5_gene_groups_axes.png",
       width = 7, height = 5)





### PCA with COG categories ----------------------------------------------

# read metadata
cog_cat_metadata = read_xlsx("data/cog_category_metadata.xlsx")

# reformat eggnog with present COG categories
eggnog_cog_cat = eggnog_COG %>% 
  drop_na(COG_category) %>% 
  filter(COG_category != "-") %>% 
  select(original_id, COG_category) %>% 
  separate_rows(COG_category, sep = "") %>% 
  filter(COG_category != "") 

eggnog_cog_cat %>% 
  write_csv("tables/eggnog_COG_categories.csv")

eggnog_cog_cat %>%
  left_join(  pca_embeddings ) %>% 
  left_join(gene_presence %>% 
              rename(original_id = Gene)) %>%
  arrange(desc(gene_group)) %>%
  ggplot(aes(x = pc1, y = pc2,
             color = COG_category)) +
  geom_point(alpha = 0.2) +
  # 2d density
  geom_density_2d(contour_var = "ndensity",
                  linewidth = 0.35,
                  color = 'grey40') +
  labs(title = "PCA projection of T5 embeddings",
       subtitle = "Color by COG category",
       x = glue("Principal component 1 ({per_var[1]})"), 
       y = glue("Principal component 2 ({per_var[2]})")) +
  facet_wrap(~COG_category)

ggsave("exploration/t5_embeddings/pca_projection_T5_COG_categories.pdf",
       width = 15, height = 12)


eggnog_cog_cat %>%
  left_join(  pca_embeddings ) %>% 
  left_join(gene_presence %>% 
              rename(original_id = Gene)) %>%
  arrange(desc(gene_group)) %>%
  drop_na(gene_group) %>%
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud"))) %>%
  # filter(gene_group == "Cloud") %>%
  ggplot(aes(x = pc1, y = pc2,
             color = gene_group)) +
  geom_point(alpha = 0.2) +
  labs(title = "PCA projection of T5 embeddings",
       subtitle = "Color by gene group",
       x = glue("Principal component 1 ({per_var[1]})"), 
       y = glue("Principal component 2 ({per_var[2]})")) +
  # 2d density
  geom_density_2d(contour_var = "ndensity",
                  linewidth = 0.35,
                  alpha = 0.3,
                  color = 'grey40') +
  scale_color_manual(values = gene_class_colors) +
  ylim(-0.8,0.7) +
  xlim(-0.7, 2)



# create a single plot per COG_category

COG_cats = eggnog_cog_cat %>% distinct(COG_category) %>% 
  arrange(COG_category) %>% pull(COG_category)
 

marginal_cog_plotter <- function(
  point_alpha = 0.3,  
  cog = "C"
    ) {
  p = eggnog_cog_cat %>%
    left_join(  pca_embeddings ) %>% 
    left_join(gene_presence %>% 
                rename(original_id = Gene)) %>%
    arrange(desc(gene_group)) %>%
    drop_na(gene_group) %>%
    mutate(gene_group = factor(gene_group, 
                               levels = c("Core", "Shell", "Cloud"))) %>%
    filter(COG_category == cog) %>%
    ggplot(aes(x = pc1, y = pc2,
               color = gene_group)) +
    geom_point(alpha = point_alpha) +
    # 2d density
    geom_density_2d(contour_var = "ndensity",
                    linewidth = 0.35,
                    alpha = 0.6,
                    aes(color = gene_group)) +
    scale_color_manual(values = gene_class_colors) +
    labs(title = glue::glue("COG category: {cog}"),
         x = glue("Principal component 1 ({per_var[1]})"), 
         y = glue("Principal component 2 ({per_var[2]})")) +
    ylim(-22,22) +
    xlim(-50, 15)
  
  
  
  p = ggMarginal(p,
                 type = "density",
                 groupColour = TRUE,
                 groupFill = TRUE) 
  
  p
}


for (cog in COG_cats) {
  print(glue("COG category: {cog}"))
  p = marginal_cog_plotter(cog = cog)
  ggsave(plot = p, 
         filename = glue::glue("exploration/t5_embeddings/PCA_COG/pca_projection_T5_COG_{cog}.pdf"),
         width = 8, height = 6)
  ggsave(plot = p,
         filename = glue::glue("exploration/t5_embeddings/PCA_COG/pca_projection_T5_COG_{cog}.png"),
         width = 8, height = 6)
}




## PC1 and PC2 analysis ---------------------------------------------------

pca_eggnog_long = eggnog_cog_cat %>%
  left_join(  pca_embeddings ) %>% 
  pivot_longer(c(pc1, pc2), 
               names_to = "PC", values_to = "value") %>% 
  mutate(PC = str_replace_all(PC, "pc1", "PC1")) %>% 
  mutate(PC = str_replace_all(PC, "pc2", "PC2")) 



pca_eggnog_long %>% 
  # filter(PC == "PC1") %>% 
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud"))) %>%
  left_join(cog_cat_metadata) %>% 
  filter(!(COG_category %in% c("B", "W", "Y"))) %>% 
  ggplot(aes(y = value, x = COG_category, fill = gene_group)) +
  geom_boxplot() +
  facet_wrap(vars(PC, Major_category), scales = "free", ncol = 4) +
  scale_fill_manual(values = gene_class_colors) +
  labs(y = "Principal component value",
       x = "COG category",
       fill = "Gene group") 


# stats
pca_cog_stats = pca_eggnog_long %>% 
  left_join(cog_cat_metadata) %>%
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud"))) %>% 
  filter(!(COG_category %in% c("B", "W", "Y"))) %>% 

  group_by(PC, COG_category, Major_category) %>%
  rstatix::wilcox_test(value ~ gene_group, 
                       detailed = TRUE,
                       p.adjust.method = "BH") %>% 
  add_xy_position(x = "COG_category", 
                  scales = "free", 
                  dodge = 1.2)

# arrange stats by Major_category and COG_category and PC
pca_cog_stats = pca_cog_stats %>% 
  arrange(Major_category, COG_category, PC)

pca_cog_stats %>% write_csv("tables/pca_T5_cog_stats.csv")

# fix coordinates for stats
pca_cog_stats = pca_cog_stats %>% 
  group_by(PC, Major_category, group1, group2) %>% 
  mutate(x = 1:n(),
         int = floor(xmin),
         xmin = xmin - int,
         xmax = xmax - int,
         xmin = (x-1) + xmin + 0.1,
         xmax = (x-1) + xmax + 0.1) %>% 
  select(-int) %>% 
  ungroup

# fix y position
pca_cog_stats = pca_cog_stats %>% 
  mutate(y.position = case_when(PC == "PC1" & group1 == "Core" & group2 == "Cloud" ~ y.position + 2.4,
                                PC == "PC1" & group1 == "Shell" & group2 == "Cloud" ~ y.position + 4.9,
                                PC == "PC2" & group1 == "Core" & group2 == "Cloud" ~ y.position + 1.1,
                                PC == "PC2" & group1 == "Shell" & group2 == "Cloud" ~ y.position + 1.1,
                                TRUE ~ y.position)) 



pca_eggnog_long %>% 
  left_join(cog_cat_metadata) %>% 
  filter(!(COG_category %in% c("B", "W", "Y"))) %>% 
  arrange(COG_category) %>% 
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud"))) %>% 
  ggpubr::ggboxplot(x = "COG_category", y = "value",
                  color = "black", 
                  fill = "gene_group",
                  facet.by = c("PC", "Major_category"),
                  # facet in 4 columns
                  scales = "free",
                  ncol = 1,
                  short.panel.labs = FALSE,
                  palette = gene_class_colors) +
  ggpubr::stat_pvalue_manual(pca_cog_stats, 
                             label = "p.adj.signif",
                             tip.length = 0.01,
                             hide.ns = F,
                             scales = "free")

ggsave("exploration/t5_embeddings/pca_projection_T5_COG_boxplot.pdf",
       width = 12, height = 8)
ggsave("exploration/t5_embeddings/pca_projection_T5_COG_boxplot.png",
       width = 12, height = 8)



### barplot of n elements per COG category, major category, gene group and PC

pca_eggnog_long %>% 
  drop_na(gene_group) %>% 
  left_join(cog_cat_metadata) %>% 
  filter(!(COG_category %in% c("B", "W", "Y"))) %>% 
  group_by(PC, COG_category, Major_category, gene_group) %>% 
  count() %>% 
  group_by(PC, COG_category, Major_category) %>%
  mutate(total_n = sum(n)) %>% 
  filter(PC == "PC1") %>% 
  ggplot(aes(y = fct_reorder(COG_category, total_n), 
             x = n, 
             fill = gene_group)) +
  geom_bar(stat = "identity") +
  # facet_wrap(vars(PC), scales = "free_x", ncol = 4) +
  scale_fill_manual(values = gene_class_colors) +
  labs(y = "COG category",
       x = "Number of elements",
       fill = "Gene group") +
  theme_cowplot(8) +
  background_grid()

ggsave("exploration/t5_embeddings/pca_projection_T5_COG_barplot.pdf",
       width = 3, height = 6)



### % of genes with COG category -------------------------

genes_with_cog = eggnog_cog_cat %>% 
  # filter(COG_category != "S") %>% 
  distinct(original_id) %>% 
  count() %>% pull(n)

(genes_with_cog / dim(gene_presence)[1]) * 100


# Strain T5 embedding average integration ------------------------------------

## Data preparation ------------------
### T5 embeddings =========================
  
# read reduced embeddings from T5 model
reduced_embeddings = read_csv("data/t5_embeddings/t5_embeddings/reduced_embeddings_file.csv") %>% 
  rename(id = `...1`) %>% 
  left_join(umap_proj %>% select(id, original_id)) %>% 
  select(-id, original_id, everything())

reduced_embeddings = reduced_embeddings %>% select(original_id, `0`:`1023`)


pca_embeddings_clusters = read_table("data/t5_embeddings/t5_embeddings/PCA_embeddings_clusters.tsv")

multicluster = pca_embeddings_clusters %>% 
  group_by(cluster) %>%
  dplyr::count() %>%
  filter(n > 1) %>% ungroup

dim(multicluster)[1]

missing_red_emb = pca_embeddings_clusters %>% 
  filter(cluster %in% multicluster$cluster) %>% 
  left_join(reduced_embeddings) %>% 
  drop_na() %>% 
  select(-original_id, -position) %>% 
  full_join(pca_embeddings_clusters %>% 
              filter(cluster %in% multicluster$cluster)) %>% 
  select(cluster, position, original_id, everything()) %>% 
  arrange(cluster) %>% 
  drop_na(original_id) %>%
  select(-cluster, -position) %>% 
  distinct(original_id, .keep_all = TRUE)


reduced_embeddings = reduced_embeddings %>%
  bind_rows(missing_red_emb) %>% 
  distinct(original_id, .keep_all = TRUE)

# there are some proteins without embeddings, make a list

missing_embeddings = reduced_embeddings %>% 
  filter(is.na(`1`)) %>% 
  pull(original_id)


core_genes = gene_presence %>% 
  filter(gene_prop > 0.95) %>% 
  pull(Gene)

gene_removals = unique(c(missing_embeddings, core_genes))


### Gene PA ----------------------------------------------------------
# read the gene PA matrix and transpose it
gene_pa = read_delim("data/panaroo_results/gene_presence_absence.Rtab", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

gene_pa = gene_pa %>% 
  filter(!(Gene %in% gene_removals)) 

gc()

reduced_embeddings = reduced_embeddings %>% 
  filter(!(original_id %in% gene_removals)) 

gene_pa$Gene
# fix missmatches 
diff1 = setdiff(reduced_embeddings$original_id, gene_pa$Gene)
diff2 = setdiff(gene_pa$Gene, reduced_embeddings$original_id)

gene_removals_missmatches = unique(c(diff1, diff2))

gene_pa = gene_pa %>% 
  filter(!(Gene %in% gene_removals_missmatches)) %>% 
  arrange(Gene)

reduced_embeddings = reduced_embeddings %>%
  filter(!(original_id %in% gene_removals_missmatches)) %>% 
  arrange(original_id)


## check that both gene sets are the same and in the same order
all(gene_pa$Gene == reduced_embeddings$original_id)
# all(names(gene_pa)[-1] == reduced_embeddings$original_id)

# transpose
gene_pa = gene_pa %>% 
  column_to_rownames("Gene") %>%
  t() %>%
  as_tibble() %>% 
  mutate(genome = genome_names$value,
         .before = "aCTx2")

gene_pa


## matrix multiplication ================================================

genes_embedd = reduced_embeddings$original_id
genomes_embedd = gene_pa$genome


embeddings = reduced_embeddings %>% 
  column_to_rownames("original_id") %>% 
  as.matrix

gene_pa_matrix = gene_pa %>% 
  column_to_rownames("genome") %>% 
  as.matrix

dim(gene_pa_matrix)
dim(embeddings)


# dimensions match for matrix multiplication
dim(gene_pa_matrix)[2] == dim(embeddings)[1]

# matrix multiplication
strain_embeddings = gene_pa_matrix %*% embeddings

dim(strain_embeddings)


colnames(strain_embeddings)
rownames(strain_embeddings)

strain_embeddings[1, ]



# save the strain embeddings
write_csv(strain_embeddings %>% 
            as_tibble(rownames = "genome"), 
          "data/t5_embeddings/t5_embeddings/strain_embeddings.csv")


# remove heavy objects
rm(embeddings, gene_pa_matrix, reduced_embeddings)
gc()



  
## PCA from strain embeddings --------------


pc = prcomp(strain_embeddings,
             center = TRUE,
             scale. = TRUE)

attributes(pc)

# percentage of variation 
var_exp_tbl = summary(pc)$importance %>% 
  as_tibble(rownames = "cosa") %>% 
  select(cosa, PC1:PC20)



pc_tb = pc %>% broom::tidy() %>% 
  filter(PC %in% c(1:10))

### FIX THIS

pca_strains = pc_tb$row
metadata_strains = metadata_final$assembly_id_simp

# differences
setdiff(pca_strains, metadata_strains)

pc_tb = pc_tb %>% 
  mutate(row = case_when(str_detect(row, "^(NT\\d{5})_\\d+") ~ str_sub(row, 1, 7),
                         TRUE ~ row))

pca_strains = pc_tb %>% 
  filter(PC %in% c(1,2)) %>% 
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>% 
  left_join(metadata_final %>%  select(row = assembly_id_simp, phylogroup)) %>% 
  drop_na(phylogroup) 


pca_strains %>% 
  ggplot(aes(x = PC1, y = PC2, color = phylogroup)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = phylo_colors) +
  theme_cowplot(14) +
  labs(title = "PCA projection of strain embeddings",
       subtitle = "Color by phylogroup",
       x = glue("Principal component 1 ({round(var_exp_tbl[2,2]*100, 2)}%)"), 
       y = glue("Principal component 2 ({round(var_exp_tbl[2,3]*100, 2)}%)"))

ggsave("exploration/strain_embeddings/pca_strain_embeddings.pdf",
       width = 10, height = 8)
ggsave("exploration/strain_embeddings/pca_strain_embeddings.png",
       width = 10, height = 8)


pca_strains %>% 
  ggplot(aes(x = PC1, y = PC2, fill = phylogroup,
             color = phylogroup)) +
  # density 2d
  geom_density_2d() +
  scale_fill_manual(values = phylo_colors) +
  scale_color_manual(values = phylo_colors) +
  theme_cowplot(14) +
  facet_wrap(~phylogroup) +
  labs(title = "PCA projection of strain embeddings",
       subtitle = "Color by phylogroup",
       x = glue("Principal component 1 ({round(var_exp_tbl[2,2]*100, 2)}%)"), 
       y = glue("Principal component 2 ({round(var_exp_tbl[2,3]*100, 2)}%)"))

ggsave("exploration/strain_embeddings/pca_strain_embeddings_density.pdf",
       width = 10, height = 8)
ggsave("exploration/strain_embeddings/pca_strain_embeddings_density.png",
       width = 10, height = 8)


### compose plot of PC1 and PC2 with phylogroups -------

pca_strains = pc_tb %>% 
  filter(PC %in% c(1,2)) %>% 
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>% 
  left_join(metadata_final %>%  select(row = assembly_id_simp, phylogroup)) %>% 
  drop_na(phylogroup) 

# Define the desired order of phylogroups
desired_phylogroup_order <- c("A", "B1", "C", 
                              "B2",
                              "D", "E", "U", "F", "G")


# Generate facet plots
phylogroups_to_plot <- desired_phylogroup_order # Use the ordered vector

plots_list <- lapply(phylogroups_to_plot, function(target_phylogroup) {
  pca_strains %>%
    mutate(highlight = ifelse(phylogroup == target_phylogroup, phylogroup, "other")) %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(data = . %>% filter(highlight == "other"), color = "grey", alpha = 0.6, size=3) +
    geom_point(data = . %>% filter(highlight != "other"), aes(color=highlight), alpha = 0.6, size = 3) +
    scale_color_manual(values = setNames(phylo_colors[target_phylogroup], target_phylogroup)) + # Here is the fix
    theme_cowplot(14) +
    labs(title = glue("Phylogroup {target_phylogroup}"),
         # subtitle = "Highlighted phylogroup",
         x = NULL, 
         y = NULL) +
    theme(legend.position = "none")
  
})

# Combine plots using patchwork and add the common axes labels
combined_plot <- wrap_plots(plots_list, ncol = 3) +
  plot_annotation(
    title = "PCA projection of strain embeddings",
    subtitle = "Highlighted phylogroup",
    theme = theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust=0.5))
  ) &
  labs(x = glue("Principal component 1"), y = glue("Principal component 2"))


print(combined_plot)

ggsave("exploration/strain_embeddings/pca_strain_embeddings_facet.pdf",
       width = 12, height = 10)
ggsave("exploration/strain_embeddings/pca_strain_embeddings_facet.png",
       width = 12, height = 10)




### boxplot of PC1 and PC2 per phylogroup ----------------------------

pca_strains_long = pc_tb %>% 
  filter(PC %in% c(1,2)) %>% 
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>% 
  left_join(metadata_final %>%  select(row = assembly_id_simp, phylogroup)) %>% 
  drop_na(phylogroup) %>% 
  pivot_longer(c(PC1, PC2), 
               names_to = "PC", values_to = "value") 


pca_strains_stats = pca_strains_long %>% 
  group_by(PC) %>% 
  rstatix::wilcox_test(value ~ phylogroup, detailed = TRUE) %>% 
  add_xy_position() 


pca_strains_long %>%
  ggpubr::ggboxplot(x = "phylogroup", y = "value",
                  color = "black", 
                  fill = "phylogroup",
                  scales = "free",
                  facet.by = "PC",
                  ncol = 1,
                  short.panel.labs = FALSE,
                  palette = phylo_colors) +
  ggpubr::stat_pvalue_manual(pca_strains_stats, 
                             label = "p.adj.signif",
                             tip.length = 0.01,
                             hide.ns = F,
                             scales = "free") +
  labs(y = "Principal Component value",
       x = NULL,
       fill = "Phylogroup") 

# calculate centroids

centroids = pca_strains_long %>%
  group_by(phylogroup, PC) %>% 
  summarise(mean = mean(value),
            std = sd(value)) %>% 
  arrange(mean)

# plot centroids with PC1 and PC2 as axes and errorbars 
centroids %>% 
  # pivot wider both mean and sd
  pivot_wider(names_from = PC, values_from = c(mean, std)) %>% 
  # plot
  ggplot(aes(x = mean_PC1, y = mean_PC2, fill = phylogroup)) +
  geom_errorbar(aes(ymin = mean_PC2 - std_PC2, ymax = mean_PC2 + std_PC2),
                width = 0.1) +
  geom_errorbarh(aes(xmin = mean_PC1 - std_PC1, xmax = mean_PC1 + std_PC1),
                 height = 0.1) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = phylo_colors) +
  labs(title = "PCA projection of strain embeddings",
       subtitle = "Centroids per phylogroup",
       x = glue("Principal component 1 ({round(var_exp_tbl[2,2]*100, 2)}%)"), 
       y = glue("Principal component 2 ({round(var_exp_tbl[2,3]*100, 2)}%)")) +
  theme_cowplot(14)

ggsave("exploration/strain_embeddings/pca_strain_embeddings_centroids.pdf",
       width = 8, height = 6)
ggsave("exploration/strain_embeddings/pca_strain_embeddings_centroids.png",
       width = 8, height = 6)


#### MANOVA from PCA --------
# using the first 10 PC because they account for 90% of variability


library(vegan)

pca_strains_10 = pc_tb %>% 
  filter(PC %in% c(1:10)) %>% 
  pivot_wider(names_from = PC, values_from = value, names_prefix = "PC") %>% 
  left_join(metadata_final %>%  select(row = assembly_id_simp, phylogroup)) %>% 
  drop_na(phylogroup) 


pca_strains_10 %>% 
  write_csv("tables/PCA_strains_embeddings_phylogroup.csv")

# 2. Perform PERMANOVA (adonis)

# First, create the dataframe with the PC columns we need
pca_matrix <- pca_strains_10 %>% select(PC1:PC2)

# adonis_result <- adonis2(pca_matrix ~ phylogroup,
#                          data = pca_strains_10,
#                          parallel = 8,
#                          permutations = 100, # permutations for p-value,
#                          method = "euclidean" #distance method
# )
# print(adonis_result)

# 3. Post-hoc Pairwise Comparisons (if overall test is significant)
library(RVAideMemoire)
pairwise_pca_strains = pairwise.perm.manova(pca_matrix,
                                            pca_strains_10$phylogroup, 
                                            test = "Pillai",
                                            nperm = 500)

pairwise_pca_strains






### Lab strains representation ===========

metadata_final %>% 
  filter(str_detect(assembly_id_simp, "NT12"))

pca_strains %>% 
  filter(str_detect(row, "NT12")) %>% 
  ggplot(aes(x = PC1, y = PC2, color = phylogroup)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = phylo_colors) +
  theme_cowplot(14) +
  labs(title = "PCA projection of strain embeddings",
       subtitle = "Color by phylogroup",
       x = glue("Principal component 1 ({round(var_exp_tbl[2,2]*100, 2)}%)"), 
       y = glue("Principal component 2 ({round(var_exp_tbl[2,3]*100, 2)}%)"))






## UMAP from strain embeddings --------------

library(umap)

n_neigh = c(10, 15, 20, 25, 30, 35, 40, 50, 60, 70)
min_dist = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)


for (n in n_neigh) {
  for (d in min_dist) {
    print(glue::glue("Running UMAP with n: {n}, d: {d}"))
    umap_strain = umap(strain_embeddings, 
                       n_neighbors = n, 
                       min_dist = d, 
                       n_components = 2)
    
    umap_strain = umap_strain$layout %>% as_tibble(rownames = "assembly_id_simp") %>% 
      rename(umap_1 = V1, umap_2 = V2) %>% 
      left_join(metadata_final %>%  select(assembly_id_simp, phylogroup)) %>% 
      drop_na(phylogroup)
    
    p = umap_strain %>%
      ggplot(aes(x = umap_1, y = umap_2, color = phylogroup)) +
      geom_point(alpha = 0.5) +
      scale_color_manual(values = phylo_colors) +
      theme_cowplot(14)
    
    ggsave(plot = p, 
           filename = glue::glue("exploration/strain_embeddings/umap_search/umap_strain_embeddings_n{n}_d{d}.pdf"),
           width = 8, height = 6)
    ggsave(plot = p, 
           filename = glue::glue("exploration/strain_embeddings/umap_search/umap_strain_embeddings_n{n}_d{d}.png"),
           width = 8, height = 6)
    print("\n")
    print("\n")
    
  }
}


umap_strain = umap(strain_embeddings, 
                   n_neighbors = 30, 
                   min_dist = 0.4, 
                   n_components = 2)

umap_strain = umap_strain$layout %>% as_tibble(rownames = "assembly_id_simp") %>% 
  rename(umap_1 = V1, umap_2 = V2) %>% 
  left_join(metadata_final %>%  select(assembly_id_simp, phylogroup)) %>% 
  drop_na(phylogroup)

umap_strain %>%
  ggplot(aes(x = umap_1, y = umap_2, color = phylogroup)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = phylo_colors) +
  theme_cowplot(14) 

ggsave("exploration/strain_embeddings/umap_strain_embeddings.pdf",
       width = 8, height = 6)
ggsave("exploration/strain_embeddings/umap_strain_embeddings.png",
       width = 8, height = 6)

# with 2d density plots
umap_strain %>%
  mutate(phylogroup = factor(phylogroup, 
                             levels = c("A", "B1",  "C",
                                        "B2",
                                        "D", "E", "U", "F", "G"))) %>%
  ggplot(aes(x = umap_1, y = umap_2, color = phylogroup)) +
  geom_point(alpha = 0.15) +
  geom_density_2d(color = 'grey10',
                  contour_var = "ndensity") +
  scale_color_manual(values = phylo_colors) +
  theme_cowplot(14) +
  facet_wrap(~phylogroup)

ggsave("exploration/strain_embeddings/umap_strain_embeddings_density.pdf",
       width = 12, height = 10)
ggsave("exploration/strain_embeddings/umap_strain_embeddings_density.png",
       width = 12, height = 10)






# genes with f(x) per method per genome -----------------------------------

gene_pa_long 


all_tidy = read_csv("tables/all_tidy_FINAL.csv") %>% 
  select(gene_family, method) %>% 
  distinct(gene_family, method) %>% 
  mutate(presence = 1)

all_tidy = all_tidy %>% pivot_wider(names_from = method, 
                                    values_from = presence,
                                    values_fill = 0) 

all_tidy = all_tidy %>%  left_join(gene_presence, by = c("gene_family" = "Gene"))

# genes present in the gene PA matrix but not in the GO tables
removals_1 = gene_pa_long %>% 
  distinct(Gene) %>% 
  pull(Gene) %>%
  setdiff(all_tidy %>% 
            select(gene_family) %>% 
            pull(gene_family))

# genes present in the table but not in the gene PA matrix
removals_2 = all_tidy %>% 
  select(gene_family) %>% 
  pull(gene_family) %>% 
  setdiff(gene_pa_long %>% 
            distinct(Gene) %>% 
            pull(Gene))

gene_removals = unique(c(removals_1, removals_2))


# filter both gene_pa_long and all_tidy
gene_pa_long_count = gene_pa_long %>% 
  filter(!(Gene %in% gene_removals)) 

all_tidy = all_tidy %>% 
  filter(!(gene_family %in% gene_removals))


# do the columns match now? 

unique(gene_pa_long_count$Gene) %>% length()
unique(all_tidy$gene_family) %>% length()

# should be 92424
sum(unique(gene_pa_long_count$Gene) %in% unique(all_tidy$gene_family))



# generate how many genes are per strain
genes_per_strain = gene_pa_long_count %>% 
  group_by(assembly_id_simp) %>% 
  count()

# histogram of genes per strain
genes_per_strain %>% 
  ggplot(aes(n)) +
  geom_histogram(binwidth = 10) +
  theme_cowplot(14) +
  labs(title = "Number of genes per strain",
       x = "Number of genes", y = "Count")

ggsave("exploration/genes_per_strain.pdf", width = 6, height = 4)

# let's calculate how many genes per strain have been annotated 
# with each method

annotation_summary_stats = gene_pa_long_count %>% 
  rename(gene_family = Gene) %>% 
  select(-presence, -phylogroup) %>% 
  left_join(all_tidy) %>% 
  group_by(assembly_id_simp) %>% 
  summarise(across(c(eggnog, interpro, proteinfer, gopredsim),
                   ~ sum(.x))) %>% 
  left_join(genes_per_strain) %>% 
  mutate(eggnog_annotated = eggnog / n,
         interpro_annotated = interpro / n,
         proteinfer_annotated = proteinfer / n,
         gopredsim_annotated = gopredsim / n) %>% 
  mutate(eggnog_unnanotated = 1 - eggnog_annotated,
         interpro_unnanotated = 1 - interpro_annotated,
         proteinfer_unnanotated = 1 - proteinfer_annotated,
         gopredsim_unnanotated = 1 - gopredsim_annotated) 


# exploration plots
# annotated genes per method per strain
library(GGally)
ggpairs(annotation_summary_stats %>% 
          select(eggnog_annotated, interpro_annotated, 
                 proteinfer_annotated, gopredsim_annotated),
        
        # set alpha to 0.3
        lower = list(continuous = wrap("points", alpha = 0.3)),
        diag = list(continuous = wrap("barDiag", alpha = 0.6)),
        # upper = list(continuous = "box", combo = "box_no_facet"),
        # axisLabels = "none",
        title = "Annotated genes per method per strain"
        )

ggsave("exploration/annotated_genes_per_method.pdf", 
       width = 12, height = 12)
ggsave("exploration/annotated_genes_per_method.png", 
       width = 12, height = 12)

# unnanotated genes 
ggpairs(annotation_summary_stats %>% 
          select(eggnog_unnanotated, interpro_unnanotated, 
                 proteinfer_unnanotated, gopredsim_unnanotated),
        
        # set alpha to 0.3
        lower = list(continuous = wrap("points", alpha = 0.3)),
        diag = list(continuous = wrap("barDiag", alpha = 0.6)),
        # upper = list(continuous = "box", combo = "box_no_facet"),
        # axisLabels = "none",
        title = "Unnanotated genes per method per strain"
)

ggsave("exploration/unnanotated_genes_per_method.pdf", 
       width = 12, height = 12)
ggsave("exploration/unnanotated_genes_per_method.png",
       width = 12, height = 12)





annotation_summary_stats_boxplots = annotation_summary_stats %>% 
  select(assembly_id_simp:n) %>% 
  mutate(eggnog_unannotated = n - eggnog,
         interpro_unannotated = n - interpro,
         proteinfer_unannotated = n - proteinfer,
         gopredsim_unannotated = n - gopredsim) %>% 
  select(-n) %>% 
  pivot_longer(cols = -assembly_id_simp) %>% 
  mutate(name = case_when(name == "eggnog" ~ "eggnog_annotated",
                   name == "interpro" ~ "interpro_annotated",
                   name == "proteinfer" ~ "proteinfer_annotated",
                   name == "gopredsim" ~ "gopredsim_annotated",
                   TRUE ~ name)) %>% 
  separate(name, into = c("method", "annotated"), sep = "_")


annotation_summary_stats_boxplots %>% 
  group_by(annotated) %>% 
  rstatix::t_test(value ~ method, detailed = TRUE, 
                  p.adjust.method = 'fdr') %>% 
  write_csv("tables/annotation_summary_stats_boxplots_ttest.csv")

# boxplots
annotation_summary_stats_boxplots %>% 
  ggplot(aes(x = method, y = value, fill = annotated)) +
  geom_boxplot(outlier.shape = NULL) +
  scale_fill_manual(values = c("grey", "black")) +
  theme_cowplot(14) +
  labs(
       x = "Method", 
       y = "Number of genes per strain") +
  facet_wrap(~annotated)
  

annotation_summary_stats_boxplots %>% 
  mutate(method = str_replace(method, "interpro", "Interproscan"),
         method = str_replace(method, "eggnog", "eggNOG"),
         method = str_replace(method, "proteinfer", "Proteinfer"),
         method = str_replace(method, "gopredsim", "goPredSim")) %>%
  mutate(method = factor(method, 
                         levels = c("Interproscan",
                                    "eggNOG",
                                    "Proteinfer",
                                    "goPredSim"))) %>% 
  ggplot(aes(x = annotated, y = value, fill = method)) +
  geom_boxplot(outlier.shape = NULL) +
  scale_fill_manual(values = colors_methods) +
  theme_cowplot(14) +
  labs(
       x = "Method", 
       y = "Number of genes per strain") +
  facet_wrap(~method, ncol = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("exploration/annotated_genes_per_method_boxplots.pdf", 
       width = 10, height = 7)
ggsave("exploration/annotated_genes_per_method_boxplots.png",
       width = 10, height = 7)


## calculate with core, shell and cloud genomes

annotation_summary_stats_group = gene_pa_long_count %>% 
  rename(gene_family = Gene) %>% 
  select(-presence, -phylogroup) %>% 
  left_join(all_tidy) %>% 
  select(-gene_prop) %>% 
  group_by(assembly_id_simp, gene_group) %>% 
  summarise(across(c(eggnog, interpro, proteinfer, gopredsim),
                   ~ sum(.x))) %>% 
  left_join(genes_per_strain) %>%
  mutate(eggnog_annotated = eggnog / n,
         interpro_annotated = interpro / n,
         proteinfer_annotated = proteinfer / n,
         gopredsim_annotated = gopredsim / n) %>%
  mutate(eggnog_unnanotated = 1 - eggnog_annotated,
         interpro_unnanotated = 1 - interpro_annotated,
         proteinfer_unnanotated = 1 - proteinfer_annotated,
         gopredsim_unnanotated = 1 - gopredsim_annotated)


annotation_summary_stats_group_boxplots = annotation_summary_stats_group %>%
  select(assembly_id_simp:n) %>% 
  mutate(eggnog_unannotated = n - eggnog,
         interpro_unannotated = n - interpro,
         proteinfer_unannotated = n - proteinfer,
         gopredsim_unannotated = n - gopredsim) %>% 
  select(-n) %>% 
  pivot_longer(cols = -c(assembly_id_simp, gene_group)) %>% 
  mutate(name = case_when(name == "eggnog" ~ "eggnog_annotated",
                   name == "interpro" ~ "interpro_annotated",
                   name == "proteinfer" ~ "proteinfer_annotated",
                   name == "gopredsim" ~ "gopredsim_annotated",
                   TRUE ~ name)) %>% 
  separate(name, into = c("method", "annotated"))


annotation_summary_stats_group_boxplots %>%
  group_by(annotated, gene_group) %>% 
  rstatix::t_test(value ~ method, detailed = TRUE, 
                  p.adjust.method = 'fdr') %>% 
  write_csv("tables/annotation_summary_stats_group_boxplots_ttest.csv")

# boxplots
annotation_summary_stats_group_boxplots %>% 
  filter(annotated == "annotated") %>% 
  mutate(method = str_replace(method, "interpro", "Interproscan"),
         method = str_replace(method, "eggnog", "eggNOG"),
         method = str_replace(method, "proteinfer", "Proteinfer"),
         method = str_replace(method, "gopredsim", "goPredSim")) %>%
  mutate(gene_group = factor(gene_group, 
                             levels = c("Core", "Shell", "Cloud"))) %>%
  ggplot(aes(x = method, y = value, fill = method)) +
  geom_boxplot(outlier.shape = NULL) +
  scale_fill_manual(values = colors_methods) +
  # scale_fill_manual(values = c("grey", "black")) +
  theme_cowplot(14) +
  labs(
       x = "Method", 
       y = "Number of genes per strain") +
  facet_wrap(~gene_group)
