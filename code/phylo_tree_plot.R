library(tidyverse)
library(tidytree)
library(treeio)
library(seqinr)
library(ggtree)
library(openxlsx)
library(treeio)
library(tidytree)
library(here)
library(readxl)
library(colorspace)
library(phytools)


theme_set(theme_light())

# load data ---------------------------------------------------------------
# check
# this is modifies the metadata once more, so it should be OK NOW, 
# DO NOT RUN AGAIN
# metadata_final = metadata_final %>% 
#   mutate(genome_id = case_when(is.na(genome_id) ~ assembly_id_simp,
#                                       TRUE ~ genome_id)) 
# metadata_final %>% write_csv('tables/MAIN_ecoli_NCBI_metadata.csv')


metadata_final

tree = read.tree('data/core_phylo_tree/core_tree.treefile')


# # # # # # # # # # # # # # # #
# general tree ------------------------------------------------------------
# # # # # # # # # # # # # # # # 


# phylogroup info
# removing file extension that is in the tree names
tree$tip.label = str_sub(tree$tip.label, end = -5)
sp_names = c(tree$tip.label)


df = data.frame(sp_names)
colnames(df) = 'label'

# fix NT labels (they are a pain...)
nt_fix = df %>% filter(str_detect(label, "NT12")) %>% 
  # remove anything after an underscore
  mutate(label_fix = str_sub(label, end = str_locate(label, "_")[1] - 1))

df = df %>% 
  left_join(nt_fix) %>%
  mutate(label_fix = case_when(is.na(label_fix) ~ label,
                   TRUE ~ label_fix)) %>% 
  select(label = label_fix) %>% as_tibble()
  

# second fix for the NT names, should be ok now
tree$tip.label = df$label
length(tree$tip.label)

# check
metadata_final %>% 
  filter(genome_id %in% tree$tip.label) %>% dim

metadata_final %>% 
  filter(genome_id %in% df$label) %>% dim

# tree %>% as_tibble() %>% full_join(df) %>% view()


# add phylogroup info and fix NT12139 strain
df = df %>% left_join(metadata_final %>% 
                        select(genome_id, phylogroup) %>%
                        rename(label = genome_id)) 
unique(df$phylogroup)


# root tree and specify groups 
grp = list(
  'A' = df %>% filter(phylogroup == 'A') %>% pull(label),
  'B1' = df %>% filter(phylogroup == 'B1') %>% pull(label),
  'B2' = df %>% filter(phylogroup == 'B2') %>% pull(label),
  'C' = df %>% filter(phylogroup == 'C') %>% pull(label),
  'G' = df %>% filter(phylogroup == 'G') %>% pull(label),
  'F' = df %>% filter(phylogroup == 'F') %>% pull(label),
  'E' = df %>% filter(phylogroup == 'E') %>% pull(label),
  'D' = df %>% filter(phylogroup == 'D') %>% pull(label),
  'U' = df %>% filter(phylogroup == 'U') %>% pull(label)
)

tree = tree %>% groupOTU(grp, 'phylogroup')


ggtree(tree, layout = 'fan', size = 0.2)+ aes(color = phylogroup) +
  geom_tiplab(size = .5, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = phylo_colors) + 
  guides(colour = guide_legend(override.aes = list(size = 6, pch = 18)))


ggsave(file = here('exploration', 'rtree_phylogroups_distances.pdf'), 
       width = 100, 
       height = 100, units = 'mm', scale = 2, device = 'pdf')

