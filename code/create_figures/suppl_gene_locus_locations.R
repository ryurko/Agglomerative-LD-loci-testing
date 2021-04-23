# PURPOSE: Create figure with comparison of AdaPT gene-set locations

library(tidyverse)
library(cowplot)
library(broom)

# Load AdaPT model data ---------------------------------------------------

pos_esnps_model_data <-
  read_csv("data/adapt_model_data/positional_esnps/rsquared25_model_data.csv")

pos_model_data <-
  read_csv("data/adapt_model_data/positional/rsquared25_model_data.csv")

# Load the AdaPT CV results -----------------------------------------------

# Make a function to do this quicker:
read_xgb_cv_results <- function(assign_type, phenotype, rsquared) {
  readRDS(
    paste0("data/adapt_results/", assign_type, "/",
           phenotype, "_rsquared", rsquared, "_xgb_cv.rds")
  )
}

# Starting with r^2 = 0.25

# Positional + eSNPs

# ASD
pos_esnps_asd_25_xgb_cv <- read_xgb_cv_results("positional_esnps", "asd", 25)

# SCZ
pos_esnps_scz_25_xgb_cv <- read_xgb_cv_results("positional_esnps", "scz", 25)

# EA
pos_esnps_ea_25_xgb_cv <- read_xgb_cv_results("positional_esnps", "ea", 25)

# Positional

# ASD
pos_asd_25_xgb_cv <- read_xgb_cv_results("positional", "asd", 25)

# SCZ
pos_scz_25_xgb_cv <- read_xgb_cv_results("positional", "scz", 25)

# EA
pos_ea_25_xgb_cv <- read_xgb_cv_results("positional", "ea", 25)


# Helper function to return chr for each AdaPT gene-set -------------------

get_adapt_disc_chr <- function(model_data, adapt_results, target_alpha = 0.05) {
  model_data[which(adapt_results$qvals <= target_alpha),] %>%
    dplyr::select(gene_set_id, gene_chr)
}


# Create dataset and plot of locations ------------------------------------

genes_loci_disc_plot <- mutate(get_adapt_disc_chr(pos_esnps_model_data,
                                                  pos_esnps_asd_25_xgb_cv),
                               phenotype = "ASD", assign_type = "Positional + eSNPs") %>%
  bind_rows(mutate(get_adapt_disc_chr(pos_esnps_model_data,
                                      pos_esnps_scz_25_xgb_cv),
                   phenotype = "SCZ", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_disc_chr(pos_esnps_model_data,
                                      pos_esnps_ea_25_xgb_cv),
                   phenotype = "EA", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_disc_chr(pos_model_data,
                                      pos_asd_25_xgb_cv),
                   phenotype = "ASD", assign_type = "Positional"),
            mutate(get_adapt_disc_chr(pos_model_data,
                                      pos_scz_25_xgb_cv),
                   phenotype = "SCZ", assign_type = "Positional"),
            mutate(get_adapt_disc_chr(pos_model_data,
                                      pos_ea_25_xgb_cv),
                   phenotype = "EA", assign_type = "Positional")) %>%
  ggplot(aes(x = as.factor(gene_chr))) +
  geom_bar(aes(fill = assign_type, color = assign_type),
           position = "identity",
           alpha = 0.3, ) +
  ggthemes::scale_fill_colorblind() +
  ggthemes::scale_color_colorblind(guide = FALSE) +
  scale_x_discrete(breaks = seq(1:22)) +
  facet_wrap(~ phenotype, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(x = "Chromosome", y = "Number of selected genes/loci",
       fill = "Assignment type") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 24),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

save_plot("figures/suppl/f_loci_location.jpg",
          genes_loci_disc_plot, ncol = 1, nrow = 3, base_asp = 2.5)
save_plot("figures/suppl/f_loci_location.pdf",
          genes_loci_disc_plot, ncol = 1, nrow = 3, base_asp = 2.5)



# Load gene-tables --------------------------------------------------------

tidy_pos_esnps_gene_gs <-
  read_csv("data/tidy_gene_gs/positional_esnps/agglom_rsquared25_ld_loci.csv")

tidy_pos_gene_gs <-
  read_csv("data/tidy_gene_gs/positional/agglom_rsquared25_ld_loci.csv")

gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Helper function to return genes and their types

get_adapt_gene_info <- function(model_data, adapt_results, gene_gs_table,
                                gene_table = gene_type_table, target_alpha = 0.05) {
  adapt_gene_sets <- model_data[which(adapt_results$qvals <= target_alpha),] %>%
    dplyr::select(gene_set_id)
  
  gene_gs_table %>%
    filter(gene_set_id %in% adapt_gene_sets$ld_loci_id) %>%
    dplyr::select(-ld_loci_id) %>%
    distinct() %>%
    inner_join(gene_table, by = "ensembl_id")
}

adapt_results_gene_info <- mutate(get_adapt_gene_info(pos_esnps_model_data,
                                                      pos_esnps_asd_25_xgb_cv,
                                                      tidy_pos_esnps_gene_gs),
                                  phenotype = "ASD", assign_type = "Positional + eSNPs") %>%
  bind_rows(mutate(get_adapt_gene_info(pos_esnps_model_data,
                                       pos_esnps_scz_25_xgb_cv,
                                       tidy_pos_esnps_gene_gs),
                   phenotype = "SCZ", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_gene_info(pos_esnps_model_data,
                                       pos_esnps_ea_25_xgb_cv,
                                       tidy_pos_esnps_gene_gs),
                   phenotype = "EA", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_gene_info(pos_model_data,
                                       pos_asd_25_xgb_cv,
                                       tidy_pos_gene_gs),
                   phenotype = "ASD", assign_type = "Positional"),
            mutate(get_adapt_gene_info(pos_model_data,
                                       pos_scz_25_xgb_cv,
                                       tidy_pos_gene_gs),
                   phenotype = "SCZ", assign_type = "Positional"),
            mutate(get_adapt_gene_info(pos_model_data,
                                       pos_ea_25_xgb_cv,
                                       tidy_pos_gene_gs),
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
         gene_chr = factor(gene_chr, levels = as.character(1:22)))

# What are the counts by positional assignment and phenotype:
adapt_results_gene_info %>%
  group_by(phenotype, assign_type) %>%
  count()
# A tibble: 6 x 3
#     Groups:   phenotype, assign_type [6]
#     phenotype assign_type            n
#     <fct>     <chr>              <int>
#   1 ASD       Positional           338
#   2 ASD       Positional + eSNPs   513
#   3 SCZ       Positional          5206
#   4 SCZ       Positional + eSNPs  5791
#   5 EA        Positional         12614
#   6 EA        Positional + eSNPs 13347

# What's the breakdown in terms of the number of genes that were individual vs
# clustered into other loci:
adapt_results_gene_info %>%
  filter(assign_type == "Positional + eSNPs") %>%
  dplyr::left_join(tidy_pos_esnps_gene_gs, by = "ensembl_id") %>%
  group_by(phenotype, ld_loci_id) %>%
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
#   1 ASD       gene           37      37     0.0721
#   2 ASD       locus          63     476     0.928
#   3 SCZ       gene          766     766     0.132
#   4 SCZ       locus        1018    5025     0.868
#   5 EA        gene         2719    2719     0.204
#   6 EA        locus        2428   10628     0.796

# Repeat for Positional:
adapt_results_gene_info %>%
  filter(assign_type == "Positional") %>%
  dplyr::left_join(tidy_pos_gene_gs, by = "ensembl_id") %>%
  group_by(phenotype, ld_loci_id) %>%
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
#   1 ASD       gene           29      29     0.0858
#   2 ASD       locus          51     309     0.914
#   3 SCZ       gene          750     750     0.144
#   4 SCZ       locus         952    4456     0.856
#   5 EA        gene         2592    2592     0.205
#   6 EA        locus        2360   10022     0.795


# Create number of genes displayed with chr relationship ------------------

# Fit a regression model for each phenotype and positional assignment:
# # implicated genes in chr ~ 0 + # genes in chr (remove intercept)

# Create a table with the number of genes per chr for each positional assignment:
pos_esnps_chr_genes <- tidy_pos_esnps_gene_gs %>%
  separate(ld_loci_id, c("chr", "id"), remove = TRUE, sep = "_") %>%
  mutate(gene_chr = as.numeric(str_remove(chr, "chr"))) %>%
  group_by(gene_chr) %>%
  summarize(n_genes = n())
pos_chr_genes <- tidy_pos_gene_gs %>%
  separate(ld_loci_id, c("chr", "id"), remove = TRUE, sep = "_") %>%
  mutate(gene_chr = as.numeric(str_remove(chr, "chr"))) %>%
  group_by(gene_chr) %>%
  summarize(n_genes = n())

# Now summarize the number of selected genes for each chr by phenotype and
# positional assignment, joining the number of total chr to regress on:
pos_esnps_chr_results <- adapt_results_gene_info %>%
  filter(assign_type == "Positional + eSNPs") %>%
  mutate(gene_chr = as.numeric(gene_chr)) %>%
  group_by(gene_chr, phenotype, assign_type) %>%
  summarize(n_selected_genes = n()) %>%
  ungroup() %>%
  left_join(pos_esnps_chr_genes, by = "gene_chr")
pos_chr_results <- adapt_results_gene_info %>%
  filter(assign_type == "Positional") %>%
  mutate(gene_chr = as.numeric(gene_chr)) %>%
  group_by(gene_chr, phenotype, assign_type) %>%
  summarize(n_selected_genes = n()) %>%
  ungroup() %>%
  left_join(pos_chr_genes, by = "gene_chr")

# Fit these separately and get the summary statistics
n_gene_chr_lm_summary <- pos_esnps_chr_results %>%
  bind_rows(pos_chr_results) %>%
  group_by(assign_type, phenotype) %>%
  do(chr_fit = tidy(lm(n_selected_genes ~ 0 + n_genes, data = .))) %>%
  unnest(chr_fit) %>%
  # Join info used for plotting the text:
  inner_join({
    pos_esnps_chr_results %>%
      bind_rows(pos_chr_results) %>%
      group_by(phenotype) %>%
      summarize(x_coord = min(n_genes),
                y_coord = max(n_selected_genes))
  }, by = c("phenotype"))

# Plot and fit regression lines without intercepts for each:
chr_lm_fit_plot <- pos_esnps_chr_results %>%
  bind_rows(pos_chr_results) %>%
  ggplot(aes(x = n_genes, y = n_selected_genes)) +
  geom_smooth(formula = "y ~ 0 + x", method = "lm") +
  geom_text(aes(label = gene_chr)) +
  facet_grid(phenotype ~ assign_type, scales = "free_y") +
  theme_bw() +
  # Add regressions summaries to top left corner:
  geom_label(data = mutate(n_gene_chr_lm_summary,
                           pval_label = paste0("p-value = ",
                                               signif(p.value, digits = 4))),
             aes(x = x_coord * 2, y = y_coord * .95,
                 label = pval_label),
             size = 5) +
  labs(x = "Total number of genes in chromosome",
       y = "Number of implicated genes in chromosome") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 24),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))
# Save
save_plot("figures/suppl/f_chr_lm_fit.jpg",
          chr_lm_fit_plot, ncol = 2, nrow = 3)

save_plot("figures/suppl/f_chr_lm_fit.pdf",
          chr_lm_fit_plot, ncol = 2, nrow = 3)

gene_type_display <- adapt_results_gene_info %>%
  ggplot(aes(x = as.factor(gene_chr))) +
  geom_bar(aes(fill = gene_biotype)) +
  #ggthemes::scale_fill_colorblind() +
  ggsci::scale_fill_npg() +
  ggsci::scale_color_npg() +
  #ggthemes::scale_color_colorblind(guide = FALSE) +
  scale_x_discrete(breaks = seq(1:22)) +
  # Add p-value labels of regression fit:
  geom_label(data = mutate(n_gene_chr_lm_summary,
                           pval_label = paste0("p-value = ",
                                               signif(p.value, digits = 4))),
             aes(x = 18, y = y_coord * .95, label = pval_label),
             size = 5) +
  facet_grid(phenotype ~ assign_type, scales = "free_y") +
  theme_bw() +
  labs(x = "Chromosome", y = "Number of implicated genes",
       fill = "Type of gene") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 24),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

# Save
save_plot("figures/suppl/f_gene_type_locations.jpg",
          gene_type_display, ncol = 2, nrow = 3)

save_plot("figures/suppl/f_gene_type_locations.pdf",
          gene_type_display, ncol = 2, nrow = 3)

