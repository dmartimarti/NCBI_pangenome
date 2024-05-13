
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



# filter metadata -----------------------------------------------------------

# only need to do this once

metadata = read_xlsx('ecoli_ncbi_metadata.xlsx') 

metadata = metadata %>% 
  filter(assembly_level %in% c("Complete Genome",
                                 "Chromosome",
                                 "Scaffold")) %>% 
  filter(checkM_marker_set == 'Escherichia coli')


write.xlsx(metadata, "ecoli_ncbi_metadata.xlsx",
           overwrite = TRUE)


## read metadata -----------------------------------------------------------
## START FROM HERE 

metadata = read_xlsx('ecoli_ncbi_metadata.xlsx') 

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

# plot distribution of scaffold_n50
metadata %>% 
  ggplot(aes(scaffold_n50)) +
  geom_histogram(binwidth = 10000) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Scaffold N50", y = "Number of genomes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


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

metadata_filt %>% 
  # plot distribution of genes 
  ggplot(aes(genes)) +
  geom_histogram() +
  # scale_x_continuous(labels = scales::comma) +
  labs(x = "Genes", y = "Number of genomes")
# filter genomes with genes > 7000 and genes < 2000

# filter genomes with seq_length > 7000000
metadata_filt = metadata_filt %>% 
  filter(seq_length < 7000000)

  
write_csv(metadata_filt, "ecoli_ncbi_metadata_final.csv")


## sequencing technology ---------------------------------------------------
metadata_filt %>% 
  separate_rows(sequencing_tech, sep = "; ") %>% 
  drop_na(sequencing_tech) %>%
  distinct(sequencing_tech) %>%
  arrange(sequencing_tech) %>% view
  
  
  
  
