
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
library(cowplot)
library(viridis)
library(here)
library(openxlsx)

theme_set(theme_cowplot(14))



# filter metadata -----------------------------------------------------------

# only need to do this once

metadata = read_xlsx('ecoli_ncbi_metadata.xlsx') 

metadata = metadata %>% 
  filter(`Assembly Level` %in% c("Complete Genome",
                                 "Chromosome",
                                 "Scaffold")) %>% 
  filter(`CheckM marker set` == 'Escherichia coli')


write.xlsx(metadata, "ecoli_ncbi_metadata.xlsx",
           overwrite = TRUE)


# read metadata -----------------------------------------------------------
## START FROM HERE 

metadata = read_xlsx('ecoli_ncbi_metadata.xlsx') 

metadata






contig_bins = seq(0, 1000, 10)

metadata %>%  
  mutate(bins = cut(contigs, 
                    breaks = contig_bins)) %>% 
  ggplot(aes(contigs, fill = bins)) +
  geom_histogram(binwidth = 10, show.legend = F)  +
  scale_fill_viridis(discrete = T)

# empty values from checkM 

metadata %>% 
  mutate(checkM_completeness_class = 
           case_when(is.na(checkM_completeness) ~ 'NA',
                     checkM_completeness < 90 ~ "< 90%",
                     checkM_completeness >= 90 ~ ">= 90%")) %>% 
  count(checkM_completeness_class)

metadata %>% 
  mutate(checkM_contamination_class = 
           case_when(is.na(checkM_contamination) ~ 'NA',
                     checkM_contamination > 5 ~ "> 5%",
                     checkM_contamination <= 5 ~ "<= 5%")) %>% 
  count(checkM_contamination_class) 


# get bioprojects from empty assemblies

metadata %>% 
  filter(is.na(assembly)) %>% 
  drop_na(bioproject) %>% 
  select(bioproject) %>% 
  distinct(bioproject) %>% 
  write_delim(here('data', "bioproject_ids.txt"))


contig_bins = seq(0, 250, 1)
metadata %>%  
  filter(checkM_completeness > 90,
         checkM_contamination < 5,
         contigs < 250) %>% 
  mutate(bins = cut(contigs, 
                    breaks = contig_bins,
                    .before = "genome_name")) %>% 
  ggplot(aes(contigs, fill = bins)) +
  geom_histogram(binwidth = 10, show.legend = F)  +
  scale_fill_viridis(discrete = T)



# filter dataset ----------------------------------------------------------

metadata %>%
  drop_na(assembly) %>% 
  mutate(checkM_completeness_class = 
           case_when(is.na(checkM_completeness) ~ 'NA',
                     checkM_completeness < 90 ~ "< 90%",
                     checkM_completeness >= 90 ~ ">= 90%")) %>% 
  count(checkM_completeness_class) %>% print

metadata %>% 
  drop_na(assembly) %>% 
  mutate(checkM_contamination_class = 
           case_when(is.na(checkM_contamination) ~ 'NA',
                     checkM_contamination > 5 ~ "> 5%",
                     checkM_contamination <= 5 ~ "<= 5%")) %>% 
  count(checkM_contamination_class) %>% print  


metadata %>% 
  drop_na(assembly) %>% 
  filter(contigs < 100) %>% 
  view




