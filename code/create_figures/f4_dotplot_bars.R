# PURPOSE: Create the figure with comparison of the AdaPT results to BH and
#          intercept-only for both positional types and all three phenotypes

library(tidyverse)
library(cowplot)

# This is just for results with LD threshold of 0.25

# Load the model data -----------------------------------------------------

pos_esnps_model_data <-
  read_csv("data/adapt_model_data/positional_esnps/rsquared25_model_data.csv")

pos_model_data <-
  read_csv("data/adapt_model_data/positional/rsquared25_model_data.csv")

# Load the intercept-only AdaPT results -----------------------------------

# Make a function to do this quicker:
read_int_only_results <- function(assign_type, phenotype, rsquared) {
  readRDS(
    paste0("data/adapt_results/", assign_type, "/",
           phenotype, "/rsquared", rsquared, "_int_only.rds")
  )
}

# Starting with r^2 = 0.25

# Positional + eSNPs

# ASD
pos_esnps_asd_25_int_only <- read_int_only_results("positional_esnps", "asd", 25)

# SCZ
pos_esnps_scz_25_int_only <- read_int_only_results("positional_esnps", "scz", 25)

# EA
pos_esnps_ea_25_int_only <- read_int_only_results("positional_esnps", "ea", 25)

# Positional

# ASD
pos_asd_25_int_only <- read_int_only_results("positional", "asd", 25)

# SCZ
pos_scz_25_int_only <- read_int_only_results("positional", "scz", 25)

# EA
pos_ea_25_int_only <- read_int_only_results("positional", "ea", 25)


# Load the AdaPT CV results -----------------------------------------------

# Make a function to do this quicker:
read_xgb_cv_results <- function(assign_type, phenotype, rsquared) {
  readRDS(
    paste0("data/adapt_results/", assign_type, "/",
           phenotype, "/rsquared", rsquared, "_xgb_cv.rds")
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


# Create dot-plots comparing number of discoveries by each ----------------

# Make a function that returns the number of BH discoveries for each vector:
return_bh_n_disc <- function(model_data_list, pval_type, target_alpha = 0.05) {
  map_dbl(model_data_list,
          function(model_data) {
            sum(p.adjust(model_data[[pval_type]],
                         method = "BH") <= target_alpha)
          })
}

return_adapt_n_disc <- function(adapt_model_list, target_alpha = 0.05) {
  map_dbl(adapt_model_list,
          function(adapt_model) {
            sum(adapt_model$qvals <= target_alpha)
          })
}

# Create the dotplot comparing the number of discoveries

adapt_dotplot <- bind_rows({
  tibble(type = c("Positional + eSNPs", "Positional"),
         BH = return_bh_n_disc(list(pos_esnps_model_data,
                                    pos_model_data),
                               pval_type = "asd_vegas_pval"),
         `AdaPT: intercept-only` =
           return_adapt_n_disc(
             list(pos_esnps_asd_25_int_only,
                  pos_asd_25_int_only)),
         `AdaPT: XGBoost` =
           return_adapt_n_disc(
             list(pos_esnps_asd_25_xgb_cv,
                  pos_asd_25_xgb_cv))) %>%
    pivot_longer(BH:`AdaPT: XGBoost`,
                 names_to = "method",
                 values_to = "n_disc") %>%
    mutate(method = fct_relevel(method, "BH"),
           type = fct_relevel(type, "Positional + eSNPs", "Positional"),
           phenotype = "ASD")
}, {
  tibble(type = c("Positional + eSNPs", "Positional"),
         BH = return_bh_n_disc(list(pos_esnps_model_data,
                                    pos_model_data),
                               pval_type = "scz_vegas_pval"),
         `AdaPT: intercept-only` =
           return_adapt_n_disc(
             list(pos_esnps_scz_25_int_only,
                  pos_scz_25_int_only)),
         `AdaPT: XGBoost` =
           return_adapt_n_disc(
             list(pos_esnps_scz_25_xgb_cv,
                  pos_scz_25_xgb_cv))) %>%
    pivot_longer(BH:`AdaPT: XGBoost`,
                 names_to = "method",
                 values_to = "n_disc") %>%
    mutate(method = fct_relevel(method, "BH"),
           type = fct_relevel(type, "Positional + eSNPs", "Positional"),
           phenotype = "SCZ")
},
{
  tibble(type = c("Positional + eSNPs", "Positional"),
         BH = return_bh_n_disc(list(pos_esnps_model_data,
                                    pos_model_data),
                               pval_type = "ea_vegas_pval"),
         `AdaPT: intercept-only` =
           return_adapt_n_disc(
             list(pos_esnps_ea_25_int_only,
                  pos_ea_25_int_only)),
         `AdaPT: XGBoost` =
           return_adapt_n_disc(
             list(pos_esnps_ea_25_xgb_cv,
                  pos_ea_25_xgb_cv))) %>%
    pivot_longer(BH:`AdaPT: XGBoost`,
                 names_to = "method",
                 values_to = "n_disc") %>%
    mutate(method = fct_relevel(method, "BH"),
           type = fct_relevel(type, "Positional + eSNPs", "Positional"),
           phenotype = "EA")
}) %>%
  mutate(phenotype = fct_relevel(phenotype, "ASD", "SCZ", "EA"),
         type = fct_rev(type)) %>%
  ggplot(aes(x = type, y = n_disc, fill = method)) +
  geom_bar(width = .2, stat = "identity",
           position = position_dodge(width = .75)) +
  geom_point(position = position_dodge(width = .75),
             size = 10, shape = 21, color = "white") +
  theme_bw() +
  ggthemes::scale_fill_colorblind() +
  ggthemes::scale_color_colorblind(guide = FALSE) +
  labs(x = "SNP-gene assignment", y = "Number of detected LD loci",
       fill = "Method") +
  coord_flip() +
  facet_wrap(~ phenotype, ncol = 3, scales = "free_x") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 26),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        plot.margin = margin(2, 2, 2, 2, "cm"))


# Next make the overlap charts -------------------------------------

# Load gene-tables
tidy_pos_esnps_gene_gs <-
  read_csv("data/tidy_gene_locus/positional_esnps/agglom_rsquared25_ld_loci.csv")

tidy_pos_gene_gs <-
  read_csv("data/tidy_gene_locus/positional/agglom_rsquared25_ld_loci.csv")

gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

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

# Generate the info about genes for each set of results
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

# Create the stack bar chart version
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
            size = 3) +
  labs(x = "Phenotype", y = "Proportion of genes") +
  scale_fill_manual(values = c("darkred", "darkorange", "darkblue")) +
  #ggsci::scale_fill_npg() +
  #ggthemes::scale_fill_colorblind() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12))

# Create an arrangement ---------------------------------------------------

# Make a modified version of the disc dot plot with two column and legend
# in the fourth spot:

adapt_dotplot_new <- adapt_dotplot +
  ggthemes::scale_fill_colorblind(guide = guide_legend(reverse = TRUE)) +
  facet_wrap(~ phenotype, ncol = 2, scales = "free_x") +
  theme(legend.position = c(.8, .25),
        legend.direction = "vertical",
        legend.key = element_rect(size = 5),
        legend.key.size = unit(3, 'lines'),
        panel.spacing = unit(2, "lines"),
        strip.text = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(y = "Number of detected genes/loci")

# Combine with the overlap bars:
dotplot_with_bars <-
  plot_grid(adapt_dotplot_new, adapt_results_overlap_bars,
            rel_widths = c(2, 1), #rel_heights = c(1.5, 1),
            labels = c("A", "B"), hjust = -.2,
            label_size = 24, label_fontface = "plain")
# Try saving this
save_plot("figures/main/f_dotplot_bars.jpg",
          dotplot_with_bars, ncol = 3, nrow = 2)
save_plot("figures/main/f_dotplot_bars.pdf",
          dotplot_with_bars, ncol = 3, nrow = 2)




