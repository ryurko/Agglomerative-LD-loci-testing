# PURPOSE: Breakdown of discoveries for higher r^2 thresholds based on status
#          of whether or not the gene was merged

library(tidyverse)


# Load model datasets -----------------------------------------------------

pos_esnps_50_model_data <-
  read_csv("data/adapt_model_data/positional_esnps/rsquared50_model_data.csv")

pos_50_model_data <-
  read_csv("data/adapt_model_data/positional/rsquared50_model_data.csv")

pos_esnps_75_model_data <-
  read_csv("data/adapt_model_data/positional_esnps/rsquared75_model_data.csv")

pos_75_model_data <-
  read_csv("data/adapt_model_data/positional/rsquared75_model_data.csv")

# Load the AdaPT CV results -----------------------------------------------

# Make a function to do this quicker:
read_xgb_cv_results <- function(assign_type, phenotype, rsquared) {
  readRDS(
    paste0("data/adapt_results/", assign_type, "/",
           phenotype, "/rsquared", rsquared, "_xgb_cv.rds")
  )
}

# Starting with r^2 = 0.50

# Positional + eSNPs

# ASD
pos_esnps_asd_50_xgb_cv <- read_xgb_cv_results("positional_esnps", "asd", 50)

# SCZ
pos_esnps_scz_50_xgb_cv <- read_xgb_cv_results("positional_esnps", "scz", 50)

# EA
pos_esnps_ea_50_xgb_cv <- read_xgb_cv_results("positional_esnps", "ea", 50)

# Positional

# ASD
pos_asd_50_xgb_cv <- read_xgb_cv_results("positional", "asd", 50)

# SCZ
pos_scz_50_xgb_cv <- read_xgb_cv_results("positional", "scz", 50)

# EA
pos_ea_50_xgb_cv <- read_xgb_cv_results("positional", "ea", 50)

# Next for r^2 = 0.75

# Positional + eSNPs

# ASD
pos_esnps_asd_75_xgb_cv <- read_xgb_cv_results("positional_esnps", "asd", 75)

# SCZ
pos_esnps_scz_75_xgb_cv <- read_xgb_cv_results("positional_esnps", "scz", 75)

# EA
pos_esnps_ea_75_xgb_cv <- read_xgb_cv_results("positional_esnps", "ea", 75)

# Positional

# ASD
pos_asd_75_xgb_cv <- read_xgb_cv_results("positional", "asd", 75)

# SCZ
pos_scz_75_xgb_cv <- read_xgb_cv_results("positional", "scz", 75)

# EA
pos_ea_75_xgb_cv <- read_xgb_cv_results("positional", "ea", 75)


# Helper function to return chr for each AdaPT gene-set -------------------

get_adapt_disc_chr <- function(model_data, adapt_results, target_alpha = 0.05) {
  model_data[which(adapt_results$qvals <= target_alpha),] %>%
    dplyr::select(ld_locus_id, chr.x) %>%
    dplyr::rename(chr = chr.x)
}


# Load gene-tables --------------------------------------------------------

tidy_pos_esnps_50_gene_gs <-
  read_csv("data/tidy_gene_locus/positional_esnps/agglom_rsquared50_ld_loci.csv")

tidy_pos_50_gene_gs <-
  read_csv("data/tidy_gene_locus/positional/agglom_rsquared50_ld_loci.csv")

tidy_pos_esnps_75_gene_gs <-
  read_csv("data/tidy_gene_locus/positional_esnps/agglom_rsquared75_ld_loci.csv")

tidy_pos_75_gene_gs <-
  read_csv("data/tidy_gene_locus/positional/agglom_rsquared75_ld_loci.csv")

gene_type_table <-
  read_csv("clean_data/gencode/gencode_v21_table.csv")

# Helper function to return genes and their types

get_adapt_gene_info <- function(model_data, adapt_results, gene_gs_table,
                                gene_table = gene_type_table, target_alpha = 0.05) {
  adapt_gene_sets <- model_data[which(adapt_results$qvals <= target_alpha),] %>%
    dplyr::select(ld_locus_id)

  gene_gs_table %>%
    filter(ld_locus_id %in% adapt_gene_sets$ld_locus_id) %>%
    dplyr::select(-ld_locus_id) %>%
    distinct() %>%
    inner_join(gene_table, by = "ensembl_id")
}

adapt_results_50_gene_info <- mutate(get_adapt_gene_info(pos_esnps_50_model_data,
                                                      pos_esnps_asd_50_xgb_cv,
                                                      tidy_pos_esnps_50_gene_gs),
                                  phenotype = "ASD", assign_type = "Positional + eSNPs") %>%
  bind_rows(mutate(get_adapt_gene_info(pos_esnps_50_model_data,
                                       pos_esnps_scz_50_xgb_cv,
                                       tidy_pos_esnps_50_gene_gs),
                   phenotype = "SCZ", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_gene_info(pos_esnps_50_model_data,
                                       pos_esnps_ea_50_xgb_cv,
                                       tidy_pos_esnps_50_gene_gs),
                   phenotype = "EA", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_gene_info(pos_50_model_data,
                                       pos_asd_50_xgb_cv,
                                       tidy_pos_50_gene_gs),
                   phenotype = "ASD", assign_type = "Positional"),
            mutate(get_adapt_gene_info(pos_50_model_data,
                                       pos_scz_50_xgb_cv,
                                       tidy_pos_50_gene_gs),
                   phenotype = "SCZ", assign_type = "Positional"),
            mutate(get_adapt_gene_info(pos_50_model_data,
                                       pos_ea_50_xgb_cv,
                                       tidy_pos_50_gene_gs),
                   phenotype = "EA", assign_type = "Positional")) %>%
  mutate(gene_biotype = fct_relevel(gene_biotype,
                                    "protein_coding", "antisense",
                                    "lnc_rna", "other"),
         gene_biotype = fct_recode(gene_biotype,
                                   `Protein coding` = "protein_coding",
                                   `Antisense` = "antisense",
                                   `lncRNA` = "lnc_rna",
                                   `Other` = "other"),
         phenotype = fct_relevel(phenotype, "ASD", "SCZ", "EA"),
         chr = factor(chr, levels = as.character(1:22)))

# What's the breakdown in terms of the number of genes that were individual vs
# clustered into other loci:
adapt_results_50_gene_info %>%
  filter(assign_type == "Positional + eSNPs") %>%
  dplyr::left_join(tidy_pos_esnps_50_gene_gs, by = "ensembl_id") %>%
  group_by(phenotype, ld_locus_id) %>%
  summarize(n_genes = n()) %>%
  ungroup() %>%
  mutate(locus_type = ifelse(n_genes > 1, "locus", "gene")) %>%
  group_by(phenotype, locus_type) %>%
  summarize(n_loci = n(), n_genes = sum(n_genes)) %>%
  mutate(prop_genes = n_genes / sum(n_genes))
# phenotype locus_type n_loci n_genes prop_genes
# <fct>     <chr>       <int>   <int>      <dbl>
#   1 ASD       gene          112     112      0.289
#   2 ASD       locus          38     275      0.711
#   3 SCZ       gene         2306    2306      0.506
#   4 SCZ       locus         499    2250      0.494
#   5 EA        gene         7054    7054      0.611
#   6 EA        locus        1098    4484      0.389

# Repeat for Positional:
adapt_results_50_gene_info %>%
  filter(assign_type == "Positional") %>%
  dplyr::left_join(tidy_pos_50_gene_gs, by = "ensembl_id") %>%
  group_by(phenotype, ld_locus_id) %>%
  summarize(n_genes = n()) %>%
  ungroup() %>%
  mutate(locus_type = ifelse(n_genes > 1, "locus", "gene")) %>%
  group_by(phenotype, locus_type) %>%
  summarize(n_loci = n(), n_genes = sum(n_genes)) %>%
  mutate(prop_genes = n_genes / sum(n_genes))
# # A tibble: 6 x 5
# # Groups:   phenotype [3]
#     phenotype locus_type n_loci n_genes prop_genes
#     <fct>     <chr>       <int>   <int>      <dbl>
#   1 ASD       gene           98      98      0.316
#   2 ASD       locus          34     212      0.684
#   3 SCZ       gene         2318    2318      0.546
#   4 SCZ       locus         453    1929      0.454
#   5 EA        gene         6882    6882      0.632
#   6 EA        locus        1025    4009      0.368


adapt_results_75_gene_info <- mutate(get_adapt_gene_info(pos_esnps_75_model_data,
                                                         pos_esnps_asd_75_xgb_cv,
                                                         tidy_pos_esnps_75_gene_gs),
                                     phenotype = "ASD", assign_type = "Positional + eSNPs") %>%
  bind_rows(mutate(get_adapt_gene_info(pos_esnps_75_model_data,
                                       pos_esnps_scz_75_xgb_cv,
                                       tidy_pos_esnps_75_gene_gs),
                   phenotype = "SCZ", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_gene_info(pos_esnps_75_model_data,
                                       pos_esnps_ea_75_xgb_cv,
                                       tidy_pos_esnps_75_gene_gs),
                   phenotype = "EA", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_gene_info(pos_75_model_data,
                                       pos_asd_75_xgb_cv,
                                       tidy_pos_75_gene_gs),
                   phenotype = "ASD", assign_type = "Positional"),
            mutate(get_adapt_gene_info(pos_75_model_data,
                                       pos_scz_75_xgb_cv,
                                       tidy_pos_75_gene_gs),
                   phenotype = "SCZ", assign_type = "Positional"),
            mutate(get_adapt_gene_info(pos_75_model_data,
                                       pos_ea_75_xgb_cv,
                                       tidy_pos_75_gene_gs),
                   phenotype = "EA", assign_type = "Positional")) %>%
  mutate(gene_biotype = fct_relevel(gene_biotype,
                                    "protein_coding", "antisense",
                                    "lnc_rna", "other"),
         gene_biotype = fct_recode(gene_biotype,
                                   `Protein coding` = "protein_coding",
                                   `Antisense` = "antisense",
                                   `lncRNA` = "lnc_rna",
                                   `Other` = "other"),
         phenotype = fct_relevel(phenotype, "ASD", "SCZ", "EA"),
         chr = factor(chr, levels = as.character(1:22)))

adapt_results_75_gene_info %>%
  filter(assign_type == "Positional + eSNPs") %>%
  dplyr::left_join(tidy_pos_esnps_75_gene_gs, by = "ensembl_id") %>%
  group_by(phenotype, ld_locus_id) %>%
  summarize(n_genes = n()) %>%
  ungroup() %>%
  mutate(locus_type = ifelse(n_genes > 1, "locus", "gene")) %>%
  group_by(phenotype, locus_type) %>%
  summarize(n_loci = n(), n_genes = sum(n_genes)) %>%
  mutate(prop_genes = n_genes / sum(n_genes))
# # A tibble: 6 x 5
# # Groups:   phenotype [3]
#     phenotype locus_type n_loci n_genes prop_genes
#     <fct>     <chr>       <int>   <int>      <dbl>
#   1 ASD       gene          179     179      0.601
#   2 ASD       locus          21     119      0.399
#   3 SCZ       gene         3297    3297      0.745
#   4 SCZ       locus         323    1127      0.255
#   5 EA        gene         9071    9071      0.808
#   6 EA        locus         652    2159      0.192

# Repeat for Positional:
adapt_results_75_gene_info %>%
  filter(assign_type == "Positional") %>%
  dplyr::left_join(tidy_pos_75_gene_gs, by = "ensembl_id") %>%
  group_by(phenotype, ld_locus_id) %>%
  summarize(n_genes = n()) %>%
  ungroup() %>%
  mutate(locus_type = ifelse(n_genes > 1, "locus", "gene")) %>%
  group_by(phenotype, locus_type) %>%
  summarize(n_loci = n(), n_genes = sum(n_genes)) %>%
  mutate(prop_genes = n_genes / sum(n_genes))
# # A tibble: 6 x 5
# # Groups:   phenotype [3]
#     phenotype locus_type n_loci n_genes prop_genes
#     <fct>     <chr>       <int>   <int>      <dbl>
#   1 ASD       gene          150     150      0.622
#   2 ASD       locus          22      91      0.378
#   3 SCZ       gene         3098    3098      0.765
#   4 SCZ       locus         274     952      0.235
#   5 EA        gene         8683    8683      0.819
#   6 EA        locus         577    1918      0.181


