# PURPOSE: Create figure with comparison of AdaPT LD loci locations

library(tidyverse)
library(cowplot)

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


# Load gene-tables --------------------------------------------------------

tidy_pos_esnps_gene_ld_loci <-
  read_csv("data/tidy_gene_ld_loci/positional_esnps/agglom_rsquared25_ld_loci.csv")

tidy_pos_gene_ld_loci <-
  read_csv("data/tidy_gene_ld_loci/positional/agglom_rsquared25_ld_loci.csv")

gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Helper function to return genes and their types

get_adapt_gene_info <- function(model_data, adapt_results, gene_ld_loci_table,
                                gene_table = gene_type_table, target_alpha = 0.05) {
  adapt_ld_loci <- model_data[which(adapt_results$qvals <= target_alpha),] %>%
    dplyr::select(ld_loci_id)
  
  gene_ld_loci_table %>%
    filter(ld_loci_id %in% adapt_ld_loci$ld_loci_id) %>%
    dplyr::select(-ld_loci_id) %>%
    distinct() %>%
    inner_join(gene_table, by = "ensembl_id")
}

adapt_results_gene_info <- mutate(get_adapt_gene_info(pos_esnps_model_data,
                                                      pos_esnps_asd_25_xgb_cv,
                                                      tidy_pos_esnps_gene_ld_loci),
                                  phenotype = "ASD", assign_type = "Positional + eSNPs") %>%
  bind_rows(mutate(get_adapt_gene_info(pos_esnps_model_data,
                                       pos_esnps_scz_25_xgb_cv,
                                       tidy_pos_esnps_gene_ld_loci),
                   phenotype = "SCZ", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_gene_info(pos_esnps_model_data,
                                       pos_esnps_ea_25_xgb_cv,
                                       tidy_pos_esnps_gene_ld_loci),
                   phenotype = "EA", assign_type = "Positional + eSNPs"),
            mutate(get_adapt_gene_info(pos_model_data,
                                       pos_asd_25_xgb_cv,
                                       tidy_pos_gene_ld_loci),
                   phenotype = "ASD", assign_type = "Positional"),
            mutate(get_adapt_gene_info(pos_model_data,
                                       pos_scz_25_xgb_cv,
                                       tidy_pos_gene_ld_loci),
                   phenotype = "SCZ", assign_type = "Positional"),
            mutate(get_adapt_gene_info(pos_model_data,
                                       pos_ea_25_xgb_cv,
                                       tidy_pos_gene_ld_loci),
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

gene_type_display <- adapt_results_gene_info %>%
  ggplot(aes(x = as.factor(gene_chr))) +
  geom_bar(aes(fill = gene_biotype)) +
  ggsci::scale_fill_npg() +
  ggsci::scale_color_npg() +
  scale_x_discrete(breaks = seq(1:22)) +
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

# Compare overlap of genes by the phenotype -------------------------------

# Given the AdaPT results table above, can summarize each by the assign type
# that it is in - and do so separately for each phenotype

adapt_results_overlap_bars <- adapt_results_gene_info %>%
  group_by(phenotype, ensembl_id) %>%
  summarize(overlap_type = ifelse(n() == 2, "Overlap",
                                  assign_type)) %>%
  ungroup() %>%
  group_by(phenotype, overlap_type) %>%
  summarize(n_genes = n()) %>%
  mutate(gene_prop = n_genes / sum(n_genes)) %>%
  ungroup() %>%
  mutate(overlap_type = fct_relevel(overlap_type,
                                    "Positional", "Overlap", "Positional + eSNPs")) %>%
  ggplot(aes(x = phenotype, y = gene_prop, fill = overlap_type)) +
  geom_bar(stat = "identity", alpha = .8) +
  geom_text(aes(label = paste0(n_genes, " genes", " (", round(gene_prop * 100, 2), "%)")),
            position = position_stack(vjust = 0.5), color = "white",
            size = 4.5) +
  labs(x = "Phenotype", y = "Proportion of genes") +
  scale_fill_manual(values = c("darkred", "darkblue", "darkorange")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12))

# Combine the two plots into a single display -----------------------------

combined_overlap_type_plot <-
  plot_grid(adapt_results_overlap_bars, gene_type_display,
            rel_widths = c(1, 2), label_size = 24,
            ncol = 2, labels = c("A", "B"), label_fontface = "plain")
save_plot("figures/main/f3_gene_overlap_type.jpg",
          combined_overlap_type_plot, ncol = 3, nrow = 3, base_asp = 2)
save_plot("figures/main/f3_gene_overlap_type.pdf",
          combined_overlap_type_plot, ncol = 3, nrow = 3, base_asp = 2)


# Compare proportion of gene types ----------------------------------------

adapt_gene_type_bars <- adapt_results_gene_info %>%
  ggplot(aes(x = assign_type, fill = gene_biotype)) +
  geom_bar(position = "fill") +
  gld_locici::scale_fill_npg() +
  theme_bw() +
  coord_flip() +
  facet_wrap(~ phenotype, ncol = 1) +
  labs(x = "SNP-gene assignment",
       y = "Proportion of genes", fill = "Type of gene") +
  theme(legend.position = "bottom", strip.background = element_blank(),
        axis.title.y = element_blank())

save_plot("figures/suppl/f_gene_type.jpg",
          adapt_gene_type_bars, ncol = 1, nrow = 1)
save_plot("figures/suppl/f_gene_type.pdf",
          adapt_gene_type_bars, ncol = 1, nrow = 1)

