
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

theme_set(theme_cowplot(14))




# read GO datasets --------------------------------------------------------

eggnog = read_csv("eggnog_tidy_FINAL.csv")
interpro = read_csv("interproscan_tidy_FINAL.csv")
proteinfer = read_csv("proteinfer_tidy_FINAL.csv")
gopredsim = read_csv("gopredsim_tidy_FINAL.csv")


# min vals to avoid log(0) -----------------------------------------------
egg_min = eggnog %>% 
  filter(evalue > 0) %>% 
  summarise(min_evalue = min(evalue)) %>% pull(min_evalue)

interpro_min = interpro %>%
  filter(evalue > 0) %>% 
  summarise(min_evalue = min(evalue)) %>% pull(min_evalue)

# tidy datasets -----------------------------------------------------------

## proteinfer -------------------------------------------------------------

proteinfer_tidy = proteinfer %>% 
  filter(confidence >= 0.4) %>% 
  filter(str_detect(predicted_label, "GO")) %>% 
  select(gene_family = sequence_name,
         GO = predicted_label,
         score = confidence) %>% 
  mutate(method = "proteinfer")


## eggnog -----------------------------------------------------------------

eggnog_tidy = eggnog %>% 
  filter(evalue < 1e-5) %>%
  select(gene_family = original_id,
         evalue, 
         GO = GOs) %>% 
  separate_rows(GO, sep = ",") %>% 
  mutate(method = "eggnog")


## interpro ---------------------------------------------------------------

interpro_tidy = interpro %>% 
  filter(evalue < 1e-5) %>%
  select(gene_family = original_id,
         GO = go_term,
         evalue) %>% 
  mutate(method = "interpro")


## gopredsim -------------------------------------------------------------


gopred_dist_long = gopredsim %>% 
  select(original_id, k_nn_1_distance, k_nn_2_distance, k_nn_3_distance) %>%
  pivot_longer(-original_id, names_to = "k_nn", values_to = "distance") %>% 
  mutate(kNN = case_when(
    k_nn == "k_nn_1_distance" ~ 1,
    k_nn == "k_nn_2_distance" ~ 2,
    k_nn == "k_nn_3_distance" ~ 3
  ))


gopred_annotations_long = gopredsim %>% 
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

gopred_tidy = gopred_long %>% 
  filter(distance < 1.05) %>% 
  select(gene_family = original_id,
         distance,
         GO = annotations) %>% 
  mutate(method = "gopredsim") 



# define scores --------------------------------------------------------

# no need to rededfine anything for the proteinfer a prior, the output
# probability is a good score
proteinfer_tidy

# negative exponential of the distance value
gopred_tidy = gopred_tidy %>% 
  mutate(score = exp(1) ** -distance)

interpro_tidy = interpro_tidy %>% 
  # avoid log(0)
  mutate(evalue = if_else(evalue == 0, interpro_min, evalue)) %>%
  mutate(score = -log10(evalue)) %>% 
  drop_na(score) 

eggnog_tidy = eggnog_tidy %>%
  # avoid log(0)
  mutate(evalue = if_else(evalue == 0, egg_min, evalue)) %>%
  mutate(score = -log10(evalue)) 



## plot distances, evalues and scores ---------------------------------


gopred_tidy %>% 
  pivot_longer(-c(gene_family, GO, method), names_to = "score", 
               values_to = "value") %>% 
  ggplot(aes(x = value)) +
  # density
  geom_density(fill = "grey", alpha = 0.5) +
  # rug
  # geom_rug() 
  facet_wrap(~score, scales = "free") 

ggsave("plots/gopredsim_density.png", width = 10, height = 10, dpi = 300)

interpro_tidy %>% 
  pivot_longer(-c(gene_family, GO, method), names_to = "score", 
               values_to = "value") %>% 
  ggplot(aes(x = value)) +
  # density
  geom_density(fill = "grey", alpha = 0.5) +
  # rug
  # geom_rug() 
  facet_wrap(~score, scales = "free")

ggsave("plots/interpro_density.png", width = 10, height = 10, dpi = 300)



eggnog_tidy %>% 
  pivot_longer(-c(gene_family, GO, method), names_to = "score", 
               values_to = "value") %>% 
  ggplot(aes(x = value)) +
  # density
  geom_density(fill = "grey", alpha = 0.5) +
  # rug
  # geom_rug() 
  facet_wrap(~score, scales = "free")  

ggsave("plots/eggnog_density.png", width = 10, height = 10, dpi = 300)


proteinfer_tidy %>% 
  pivot_longer(-c(gene_family, GO, method), names_to = "score", 
               values_to = "value") %>% 
  ggplot(aes(x = value)) +
  # density
  geom_density(fill = "grey", alpha = 0.5) +
  # rug
  # geom_rug() 
  facet_wrap(~score, scales = "free")

ggsave("plots/proteinfer_density.png", width = 5, height = 10 , dpi = 300)






# Rank assingment ---------------------------------------------------------

## combine all datasets --------------------------------------------------

all_tidy = bind_rows(proteinfer_tidy, eggnog_tidy, 
                     interpro_tidy, gopred_tidy) %>% 
  mutate(score = as.numeric(score)) %>% 
  select(-evalue, -distance) %>% 
  drop_na(score, GO) %>% 
  filter(score > 0) 


## rank assignment --------------------------------------------------------

all_tidy = all_tidy %>% 
  group_by(gene_family, method) %>% 
  arrange(desc(score)) %>%
  # THIS NEEDS TO BE DESCRIBED IN THE METHODS SECTIOJN
  mutate(rank = rank(-score, ties.method = "min")) %>% 
  ungroup() 

all_tidy = all_tidy %>% 
  group_by(gene_family, GO) %>%
  mutate(min_rank = min(rank),
            median_rank = median(rank),
            max_rank = max(rank),
            average_rank = mean(rank),
            sum_rank = sum(rank)) %>% 
  ungroup %>% 
  mutate(inv_sum_rank = 1 / sum_rank)

all_tidy %>% 
  write_csv("all_tidy_FINAL.csv")

# how many times a GO term is repeated per gene
all_tidy = all_tidy %>% 
  group_by(gene_family, GO) %>%
  mutate(n_GO = n()) %>% 
  ungroup

all_tidy %>% 
  ungroup %>% 
  count(n_GO)



all_ranks = all_tidy %>% 
  select(-rank) %>% 
  distinct(gene_family, GO, .keep_all = TRUE) 


all_tidy %>% 
  distinct(GO) %>% 
  write_csv("GO_terms_unique.csv")





# INFORMATION CONTENT -----------------------------------------------------

go_unique = all_tidy %>% 
  distinct(GO)

descendants = read_csv("descendants_output.csv") %>% 
  filter(Descendant_GO %in% go_unique$GO) 

ancestors = read_csv("ancestors_output.csv") %>% 
  filter(!(str_detect(Ancestor_GO, "http://purl.obolibrary.org"))) %>% 
  filter(Ancestor_GO %in% go_unique$GO)


descendants_n = descendants %>% 
  group_by(Query_GO) %>% 
  count() %>% 
  ungroup %>% 
  mutate(n = n -1) # remove the self count

ancestors_n = ancestors %>%
  group_by(Query_GO) %>% 
  count() %>% 
  ungroup 



setdiff(go_unique$GO, descendants_n$Query_GO)
setdiff(go_unique$GO, ancestors_n$Query_GO)


# plot the distribution of the number of descendants and ancestors per GO term

descendants_n %>% 
  ggplot(aes(x = n)) +
  geom_histogram(fill = "grey", alpha = 0.8) +
  labs(title = "Number of descendants per GO term",
       x = "Number of descendants",
       y = "Frequency")

ancestors_n %>%
  ggplot(aes(x = n)) +
  geom_histogram(fill = "grey", alpha = 0.8) +
  labs(title = "Number of ancestors per GO term",
       x = "Number of ancestors",
       y = "Frequency")


## calculate p(t) and IC ----------------------------

# p(t) = 1 - (descendants / ancestors + descendants)
# p(t)_inv = abs(p(t) / max(p(t)))
# IC = -log2(p(t)_inv)
# this way we get the IC of the GO term, the higher the IC, the more specific the term
pt = descendants_n %>% 
  rename(descendants = n) %>% 
  left_join(ancestors_n %>% rename(ancestor = n)) %>% 
  mutate(descendants = as.numeric(descendants),
         ancestor = as.numeric(ancestor)) %>%
  # mutate(ancestor = ancestor - 1) %>% 
  mutate(p_t =  (descendants / (ancestor + descendants))) 

# for the cases where no descendant is present, correct it by adding 1/ancestor
no_descendants = pt %>% 
  filter(descendants == 0) %>% 
  mutate(p_t = (p_t + (1/ancestor))) 



pt = pt %>% 
  filter(descendants > 0) %>%
  bind_rows(no_descendants)

# calculate the IC
pt = pt %>% 
  mutate(IC = -log2(p_t),
         p_t_norm = abs(p_t / sum(p_t)),
         IC_shannon = -(p_t_norm * log2(p_t_norm)))




pt %>% 
  ggplot(aes(x = IC)) +
  geom_histogram(fill = "grey", alpha = 1,
                 color = 'black') +
  labs(title = "Information content distribution",
       x = "Information content",
       y = "Frequency")

ggsave("plots/IC_distribution.png", width = 10, height = 7, dpi = 300)
ggsave("plots/IC_distribution.pdf", width = 10, height = 7, dpi = 300)

pt %>% 
  write_csv("IC_GOterms.csv")




