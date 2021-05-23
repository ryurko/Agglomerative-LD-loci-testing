# PURPOSE: Generate the ASD gene-set enrichment results at the gene/locus level.
#          Will generate results unadjusted and adjusted for size for each
#          phenotype and positional assignment.

library(tidyverse)
library(cowplot)
library(fgsea)

# Load the gene/locus level results ---------------------------------------

# Start with the results for r^2 = 0.25 for both positional types:
pos_esnps_pval_data <-
  map_dfr(list.files("data/locus_level_pvals/positional_esnps/rsquared25/",
                     full.names = TRUE),
          function(file_name) read_csv(file_name))

pos_only_pval_data <-
  map_dfr(list.files("data/locus_level_pvals/positional/rsquared25/",
                     full.names = TRUE),
          function(file_name) read_csv(file_name))


# Load the tidy gene to genes/loci table ----------------------------------

tidy_pos_esnps_gene_table <-
  read_csv("data/tidy_gene_locus/positional_esnps/agglom_rsquared25_ld_loci.csv")

tidy_pos_gene_table <-
  read_csv("data/tidy_gene_locus/positional/agglom_rsquared25_ld_loci.csv")

# Load gencode data to join over the start and end positions for computing gene-set length:
gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Join the start and end positions:
tidy_pos_esnps_gene_table <- tidy_pos_esnps_gene_table %>%
  dplyr::left_join(dplyr::select(gene_type_table, ensembl_id, start, end),
                   by = "ensembl_id")
tidy_pos_gene_table <- tidy_pos_gene_table %>%
  dplyr::left_join(dplyr::select(gene_type_table, ensembl_id, start, end),
                   by = "ensembl_id")

# Load curated gene lists for GSEA ----------------------------------------

curated_gene_table <-
  readxl::read_excel("data/gsea/input/curated_gene_lists_04-06-2021.xlsx")

# Join gene list info over to make genes/loci level lists -----------------

tidy_pos_esnps_gene_table <- tidy_pos_esnps_gene_table %>%
  dplyr::left_join(curated_gene_table, by = c("ensembl_id" = "ensg"))
tidy_pos_gene_table <- tidy_pos_gene_table %>%
  dplyr::left_join(curated_gene_table, by = c("ensembl_id" = "ensg"))

# How many of each type are in the different curated gene lists

# Positional + eSNPs

# Brain expressed:
length(which(tidy_pos_esnps_gene_table$brain == "TRUE"))
# [1] 16909

# Synaptic:
length(which(tidy_pos_esnps_gene_table$synaptic == "TRUE"))
# [1] 1156

# Epigenetic:
length(which(tidy_pos_esnps_gene_table$epigenetic == "TRUE"))
# [1] 722

# Transcription factor:
length(which(tidy_pos_esnps_gene_table$transcription_factor == "TRUE"))
# [1] 1716

# ASC 185:
length(which(tidy_pos_esnps_gene_table$asc185 == "TRUE"))
# [1] 185

# ASC 102:
length(which(tidy_pos_esnps_gene_table$asc102 == "TRUE"))
# [1] 102


# Positional

# Brain expressed:
length(which(tidy_pos_gene_table$brain == "TRUE"))
# [1] 16754

# Synaptic:
length(which(tidy_pos_gene_table$synaptic == "TRUE"))
# [1] 1152

# Epigenetic:
length(which(tidy_pos_gene_table$epigenetic == "TRUE"))
# [1] 715

# Transcription factor:
length(which(tidy_pos_gene_table$transcription_factor == "TRUE"))
# [1] 1701

# ASC 185:
length(which(tidy_pos_gene_table$asc185 == "TRUE"))
# [1] 185

# ASC 102:
length(which(tidy_pos_gene_table$asc102 == "TRUE"))
# [1] 102


# Summarize gene/locus level membership in the six different lists:
pos_esnps_gs_level_list_info <- tidy_pos_esnps_gene_table %>%
  mutate(brain = ifelse(is.na(brain), FALSE, as.logical(brain)),
         synaptic = ifelse(is.na(synaptic), FALSE, as.logical(synaptic)),
         epigenetic = ifelse(is.na(epigenetic), FALSE, as.logical(epigenetic)),
         transcription_factor = ifelse(is.na(transcription_factor),
                                       FALSE, as.logical(transcription_factor)),
         asc102 = ifelse(is.na(asc102), FALSE, as.logical(asc102)),
         asc185 = ifelse(is.na(asc185), FALSE, as.logical(asc185))) %>%
  group_by(ld_locus_id) %>%
  summarize(any_brain = any(brain),
            any_synaptic = any(synaptic),
            any_epigenetic = any(epigenetic),
            any_tf = any(transcription_factor),
            any_asc102 = any(asc102),
            any_asc185 = any(asc185),
            # Compute the start and end position as the min and max across genes:
            start_pos = min(start),
            end_pos = max(end)) %>%
  mutate(length_size = end_pos - start_pos)

# Number of each type post merging
length(which(pos_esnps_gs_level_list_info$any_brain))
# [1] 11558
length(which(pos_esnps_gs_level_list_info$any_synaptic))
# [1] 1114
length(which(pos_esnps_gs_level_list_info$any_epigenetic))
# [1] 664
length(which(pos_esnps_gs_level_list_info$any_tf))
# [1] 1454
length(which(pos_esnps_gs_level_list_info$any_asc185))
# [1] 183
length(which(pos_esnps_gs_level_list_info$any_asc102))
# [1] 102

# Repeat for Positional:
pos_gs_level_list_info <- tidy_pos_gene_table %>%
  mutate(brain = ifelse(is.na(brain), FALSE, as.logical(brain)),
         synaptic = ifelse(is.na(synaptic), FALSE, as.logical(synaptic)),
         epigenetic = ifelse(is.na(epigenetic), FALSE, as.logical(epigenetic)),
         transcription_factor = ifelse(is.na(transcription_factor),
                                       FALSE, as.logical(transcription_factor)),
         asc102 = ifelse(is.na(asc102), FALSE, as.logical(asc102)),
         asc185 = ifelse(is.na(asc185), FALSE, as.logical(asc185))) %>%
  group_by(ld_locus_id) %>%
  summarize(any_brain = any(brain),
            any_synaptic = any(synaptic),
            any_epigenetic = any(epigenetic),
            any_tf = any(transcription_factor),
            any_asc102 = any(asc102),
            any_asc185 = any(asc185),
            # Compute the start and end position as the min and max across genes:
            start_pos = min(start),
            end_pos = max(end)) %>%
  mutate(length_size = end_pos - start_pos)

# Number of each type post merging
length(which(pos_gs_level_list_info$any_brain))
# [1] 11878
length(which(pos_gs_level_list_info$any_synaptic))
# [1] 1109
length(which(pos_gs_level_list_info$any_epigenetic))
# [1] 654
length(which(pos_gs_level_list_info$any_tf))
# [1] 1453
length(which(pos_gs_level_list_info$any_asc185))
# [1] 183
length(which(pos_gs_level_list_info$any_asc102))
# [1] 102

# Join these indicators over to the genes/loci results table:
pos_esnps_results_data <- pos_esnps_pval_data %>%
  dplyr::left_join(pos_esnps_gs_level_list_info,
                   by = "ld_locus_id") %>%
  # Transform p-values to one-sided z stats:
  mutate(scz_z = qnorm(1 - scz_quad_pval),
         asd_z = qnorm(1 - asd_quad_pval),
         ea_z = qnorm(1 - ea_quad_pval),
         cap_scz_quad_pval = pmax(1e-15, pmin(1 - 1e-15, scz_quad_pval)),
         cap_asd_quad_pval = pmax(1e-15, pmin(1 - 1e-15, asd_quad_pval)),
         cap_ea_quad_pval = pmax(1e-15, pmin(1 - 1e-15, ea_quad_pval)),
         cap_scz_z = qnorm(1 - cap_scz_quad_pval),
         cap_asd_z = qnorm(1 - cap_asd_quad_pval),
         cap_ea_z = qnorm(1 - cap_ea_quad_pval))

pos_only_results_data <- pos_only_pval_data %>%
  dplyr::left_join(pos_gs_level_list_info,
                   by = "ld_locus_id") %>%
  # Transform p-values to one-sided z stats:
  mutate(scz_z = qnorm(1 - scz_quad_pval),
         asd_z = qnorm(1 - asd_quad_pval),
         ea_z = qnorm(1 - ea_quad_pval),
         cap_scz_quad_pval = pmax(1e-15, pmin(1 - 1e-15, scz_quad_pval)),
         cap_asd_quad_pval = pmax(1e-15, pmin(1 - 1e-15, asd_quad_pval)),
         cap_ea_quad_pval = pmax(1e-15, pmin(1 - 1e-15, ea_quad_pval)),
         cap_scz_z = qnorm(1 - cap_scz_quad_pval),
         cap_asd_z = qnorm(1 - cap_asd_quad_pval),
         cap_ea_z = qnorm(1 - cap_ea_quad_pval))


# Create figure of z-stats against size for each phenotype and ass --------

# Pipeline to make results long with each phenotype:
make_phenotype_long <- . %>%
  dplyr::select(ld_locus_id, length_size, cap_asd_z, cap_scz_z, cap_ea_z) %>%
  pivot_longer(cap_asd_z:cap_ea_z,
               names_to = "phenotype", values_to = "z_stat") %>%
  mutate(phenotype = toupper(str_remove_all(phenotype, "cap_|_z")))

# Get regression summaries for each
size_lm_summary <- pos_esnps_results_data %>%
  make_phenotype_long() %>%
  mutate(assign_type = "Positional + eSNPs") %>%
  bind_rows(pos_only_results_data %>%
              make_phenotype_long() %>%
              mutate(assign_type = "Positional")) %>%
  mutate(phenotype = fct_relevel(phenotype, "ASD", "SCZ", "EA")) %>%
  group_by(assign_type, phenotype) %>%
  do(size_fit = broom::tidy(lm(z_stat ~ log(length_size), data = .))) %>%
  unnest(size_fit)


size_lm_plot <- pos_esnps_results_data %>%
  make_phenotype_long() %>%
  mutate(assign_type = "Positional + eSNPs") %>%
  bind_rows(pos_only_results_data %>%
              make_phenotype_long() %>%
              mutate(assign_type = "Positional")) %>%
  mutate(phenotype = fct_relevel(phenotype, "ASD", "SCZ", "EA")) %>%
  ggplot(aes(x = log(length_size), y = z_stat)) +
  geom_point(alpha = 0.10, color = "darkblue") +
  geom_smooth(method = "lm", color = "red") +
  # Add regressions summaries to top left corner:
  geom_label(data = mutate(filter(size_lm_summary,
                                  term == "log(length_size)"),
                           pval_label = paste0("p-value = ",
                                               signif(p.value, digits = 4))),
             aes(x = 14, y = -3, label = pval_label),
             size = 4) +
  facet_grid(phenotype ~ assign_type, scales = "free_y") +
  theme_bw() +
  labs(y = "Unadjusted gene/locus z statistic",
       x = "log(gene/locus size)") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 24),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

# Save
save_plot("figures/suppl/f_size_lm_fit.jpg",
          size_lm_plot, ncol = 2, nrow = 3)


# Get the size adjusted z stats for each phenotype ------------------------

# First the Positional + eSNPs regressions:
pos_esnps_size_asd_lm <- lm(cap_asd_z ~ log(length_size), data = pos_esnps_results_data)
pos_esnps_size_scz_lm <- lm(cap_scz_z ~ log(length_size), data = pos_esnps_results_data)
pos_esnps_size_ea_lm <- lm(cap_ea_z ~ log(length_size), data = pos_esnps_results_data)


# Now make the adjusted version using the intercept and residual:
pos_esnps_results_data <- pos_esnps_results_data %>%
  mutate(adj_asd_z = coef(pos_esnps_size_asd_lm)[1] + residuals(pos_esnps_size_asd_lm),
         adj_scz_z = coef(pos_esnps_size_scz_lm)[1] + residuals(pos_esnps_size_scz_lm),
         adj_ea_z = coef(pos_esnps_size_ea_lm)[1] + residuals(pos_esnps_size_ea_lm))

# Repeat for Positional:
pos_only_size_asd_lm <- lm(cap_asd_z ~ log(length_size), data = pos_only_results_data)
pos_only_size_scz_lm <- lm(cap_scz_z ~ log(length_size), data = pos_only_results_data)
pos_only_size_ea_lm <- lm(cap_ea_z ~ log(length_size), data = pos_only_results_data)
pos_only_results_data <- pos_only_results_data %>%
  mutate(adj_asd_z = coef(pos_only_size_asd_lm)[1] + residuals(pos_only_size_asd_lm),
         adj_scz_z = coef(pos_only_size_scz_lm)[1] + residuals(pos_only_size_scz_lm),
         adj_ea_z = coef(pos_only_size_ea_lm)[1] + residuals(pos_only_size_ea_lm))


# Perform GSEA using unadjusted and adjusted tests ------------------------

# Make a named list of the type of gene sets for each assignment type:
pos_esnps_target_sets <-
  list("brain" = pull(filter(pos_esnps_results_data, any_brain),
                      ld_locus_id),
       "synaptic" = pull(filter(pos_esnps_results_data, any_synaptic),
                         ld_locus_id),
       "epigenetic" = pull(filter(pos_esnps_results_data, any_epigenetic),
                           ld_locus_id),
       "tf" = pull(filter(pos_esnps_results_data, any_tf),
                   ld_locus_id),
       "asc102" = pull(filter(pos_esnps_results_data, any_asc102),
                       ld_locus_id))
pos_only_target_sets <-
  list("brain" = pull(filter(pos_only_results_data, any_brain),
                      ld_locus_id),
       "synaptic" = pull(filter(pos_only_results_data, any_synaptic),
                         ld_locus_id),
       "epigenetic" = pull(filter(pos_only_results_data, any_epigenetic),
                           ld_locus_id),
       "tf" = pull(filter(pos_only_results_data, any_tf),
                   ld_locus_id),
       "asc102" = pull(filter(pos_only_results_data, any_asc102),
                       ld_locus_id))


# Start with ASD ----------------------------------------------------------

# Unadjusted results ------------------------------------------------------

# Named vector of unadjusted z stats for ASD:
pos_esnps_unadj_asd_z <- pos_esnps_results_data$cap_asd_z
names(pos_esnps_unadj_asd_z) <- pos_esnps_results_data$ld_locus_id
pos_only_unadj_asd_z <- pos_only_results_data$cap_asd_z
names(pos_only_unadj_asd_z) <- pos_only_results_data$ld_locus_id

# Perform the unadjusted GSEA for Positional + eSNPs
set.seed(1971)
pos_esnps_asd_unadj_gsea <- fgsea(pos_esnps_target_sets,
                                  pos_esnps_unadj_asd_z,
                                  # Use 10000 gene-set permutations
                                  nperm = 10000,
                                  # Positive enrichment
                                  scoreType = "pos")

# Now for Positional
set.seed(1903)
pos_only_asd_unadj_gsea <- fgsea(pos_only_target_sets,
                                 pos_only_unadj_asd_z, nperm = 10000,
                                 scoreType = "pos")

# Now create and save a plot with the -log10 of the adjusted p-values:
simplify_gsea_results <- . %>%
  as_tibble() %>%
  dplyr::select(pathway, padj) %>%
  mutate(pathway = fct_recode(pathway,
                              `Brain expressed genes` = "brain",
                              `Synaptic` = "synaptic",
                              `Epigenetic` = "epigenetic",
                              `Transcription factor genes` = "tf",
                              `Satterstrom 102 ASD genes` = "asc102"))

unadj_asd_gsea_results <- pos_esnps_asd_unadj_gsea %>%
  simplify_gsea_results() %>%
  mutate(assign_type = "Positional + eSNPs") %>%
  bind_rows(
    pos_only_asd_unadj_gsea %>%
      simplify_gsea_results() %>%
      mutate(assign_type = "Positional")
  ) %>%
  group_by(pathway) %>%
  mutate(total_logpval = sum(-log10(padj))) %>%
  ungroup() %>%
  mutate(pathway = fct_reorder(pathway, total_logpval)) %>%
  ggplot(aes(x = pathway, y = -log10(padj),
             color = assign_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = .75),
           width = 0.1, aes(fill = assign_type)) +
  geom_point(position = position_dodge(width = .75),
             size = 6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray", size = 1.5) +
  labs(y = "-log10(BH adjusted p-value)",
       x = "Gene/locus Set") +
  scale_fill_manual(values = c("darkred", "darkblue"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = c("darkred", "darkblue"), guide = FALSE) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(.85, .15),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# Adjusted results ------------------------------------------------------

# Named vector of adjusted z stats for ASD:
pos_esnps_adj_asd_z <- pos_esnps_results_data$adj_asd_z
names(pos_esnps_adj_asd_z) <- pos_esnps_results_data$ld_locus_id
pos_only_adj_asd_z <- pos_only_results_data$adj_asd_z
names(pos_only_adj_asd_z) <- pos_only_results_data$ld_locus_id

# Perform the adjusted GSEA for Positional + eSNPs
set.seed(1909)
pos_esnps_asd_adj_gsea <- fgsea(pos_esnps_target_sets,
                                pos_esnps_adj_asd_z, nperm = 10000,
                                scoreType = "pos")

# Now for Positional
set.seed(1916)
pos_only_asd_adj_gsea <- fgsea(pos_only_target_sets,
                               pos_only_adj_asd_z, nperm = 10000,
                               scoreType = "pos")

# Create and save the tables of results, by first making a helper pipeline to use
make_gsea_table <- . %>%
  as_tibble() %>%
  dplyr::select(pathway, pval, padj) %>%
  dplyr::rename(gene_list = pathway, gsea_pval = pval,
                bh_adj_pval = padj) %>%
  mutate(gene_list = fct_recode(gene_list,
                                `Brain expressed genes` = "brain",
                                `Synaptic` = "synaptic",
                                `Epigenetic` = "epigenetic",
                                `Transcription factor genes` = "tf",
                                `Satterstrom 102 ASD genes` = "asc102"))

pos_esnps_asd_adj_gsea %>%
  make_gsea_table %>%
  mutate(phenotype = "ASD", positional_type = "Positional + eSNPs") %>%
  dplyr::select(phenotype, positional_type, everything()) %>%
  write_csv("data/gsea/asd_positional_esnps_results.csv")

pos_only_asd_adj_gsea %>%
  make_gsea_table %>%
  mutate(phenotype = "ASD", positional_type = "Positional") %>%
  dplyr::select(phenotype, positional_type, everything()) %>%
  write_csv("data/gsea/asd_positional_results.csv")


# Now create and save a plot with the -log10 of the adjusted p-values:
adj_asd_gsea_results <- pos_esnps_asd_adj_gsea %>%
  simplify_gsea_results() %>%
  mutate(assign_type = "Positional + eSNPs") %>%
  bind_rows(
    pos_only_asd_adj_gsea %>%
      simplify_gsea_results() %>%
      mutate(assign_type = "Positional")
  ) %>%
  group_by(pathway) %>%
  mutate(total_logpval = sum(-log10(padj))) %>%
  ungroup() %>%
  mutate(pathway = fct_rev(fct_relevel(pathway, "Synaptic", "Epigenetic",
                               "Transcription factor genes",
                               "Satterstrom 102 ASD genes",
                               "Brain expressed genes"))) %>%
  ggplot(aes(x = pathway, y = -log10(padj),
             color = assign_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = .75),
           width = 0.1, aes(fill = assign_type)) +
  geom_point(position = position_dodge(width = .75),
             size = 6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray", size = 1.5) +
  labs(y = "-log10(BH adjusted p-value)",
       x = "Gene/locus set") +
  scale_fill_manual(values = c("darkred", "darkblue"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = c("darkred", "darkblue"), guide = FALSE) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(.85, .15),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))


# Next with SCZ ----------------------------------------------------------

# Unadjusted results ------------------------------------------------------

# Named vector of unadjusted z stats for scz:
pos_esnps_unadj_scz_z <- pos_esnps_results_data$cap_scz_z
names(pos_esnps_unadj_scz_z) <- pos_esnps_results_data$ld_locus_id
pos_only_unadj_scz_z <- pos_only_results_data$cap_scz_z
names(pos_only_unadj_scz_z) <- pos_only_results_data$ld_locus_id

# Perform the unadjusted GSEA for Positional + eSNPs
set.seed(1979)
pos_esnps_scz_unadj_gsea <- fgsea(pos_esnps_target_sets,
                                  pos_esnps_unadj_scz_z, nperm = 10000,
                                  scoreType = "pos")

# Now for Positional
set.seed(1912)
pos_only_scz_unadj_gsea <- fgsea(pos_only_target_sets,
                                 pos_only_unadj_scz_z, nperm = 10000,
                                 scoreType = "pos")

# Adjusted results ------------------------------------------------------

# Named vector of adjusted z stats for scz:
pos_esnps_adj_scz_z <- pos_esnps_results_data$adj_scz_z
names(pos_esnps_adj_scz_z) <- pos_esnps_results_data$ld_locus_id
pos_only_adj_scz_z <- pos_only_results_data$adj_scz_z
names(pos_only_adj_scz_z) <- pos_only_results_data$ld_locus_id

# Perform the adjusted GSEA for Positional + eSNPs
set.seed(1925)
pos_esnps_scz_adj_gsea <- fgsea(pos_esnps_target_sets,
                                pos_esnps_adj_scz_z, nperm = 10000,
                                scoreType = "pos")

# Now for Positional
set.seed(1918)
pos_only_scz_adj_gsea <- fgsea(pos_only_target_sets,
                               pos_only_adj_scz_z, nperm = 10000,
                               scoreType = "pos")

# Save tables of results:
pos_esnps_scz_adj_gsea %>%
  make_gsea_table %>%
  mutate(phenotype = "SCZ", positional_type = "Positional + eSNPs") %>%
  dplyr::select(phenotype, positional_type, everything()) %>%
  write_csv("data/gsea/scz_positional_esnps_results.csv")

pos_only_scz_adj_gsea %>%
  make_gsea_table %>%
  mutate(phenotype = "SCZ", positional_type = "Positional") %>%
  dplyr::select(phenotype, positional_type, everything()) %>%
  write_csv("data/gsea/scz_positional_results.csv")

# Now create a plot with the -log10 of the adjusted p-values:
adj_scz_gsea_results <- pos_esnps_scz_adj_gsea %>%
  simplify_gsea_results() %>%
  mutate(assign_type = "Positional + eSNPs") %>%
  bind_rows(
    pos_only_scz_adj_gsea %>%
      simplify_gsea_results() %>%
      mutate(assign_type = "Positional")
  ) %>%
  group_by(pathway) %>%
  mutate(total_logpval = sum(-log10(padj))) %>%
  ungroup() %>%
  mutate(pathway = fct_reorder(pathway, total_logpval)) %>%
  ggplot(aes(x = pathway, y = -log10(padj),
             color = assign_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = .75),
           width = 0.1, aes(fill = assign_type)) +
  geom_point(position = position_dodge(width = .75),
             size = 6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray", size = 1.5) +
  labs(y = "-log10(BH adjusted p-value)",
       x = "Gene/locus set") +
  scale_fill_manual(values = c("darkred", "darkblue"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = c("darkred", "darkblue"), guide = FALSE) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(.85, .15),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# And with EA ----------------------------------------------------------

# Unadjusted results ------------------------------------------------------

# Named vector of unadjusted z stats for ea:
pos_esnps_unadj_ea_z <- pos_esnps_results_data$cap_ea_z
names(pos_esnps_unadj_ea_z) <- pos_esnps_results_data$ld_locus_id
pos_only_unadj_ea_z <- pos_only_results_data$cap_ea_z
names(pos_only_unadj_ea_z) <- pos_only_results_data$ld_locus_id

# Perform the unadjusted GSEA for Positional + eSNPs
set.seed(1960)
pos_esnps_ea_unadj_gsea <- fgsea(pos_esnps_target_sets,
                                 pos_esnps_unadj_ea_z, nperm = 10000,
                                 scoreType = "pos")

# Now for Positional
set.seed(1915)
pos_only_ea_unadj_gsea <- fgsea(pos_only_target_sets,
                                pos_only_unadj_ea_z, nperm = 10000,
                                scoreType = "pos")

# Adjusted results ------------------------------------------------------

# Named vector of adjusted z stats for ea:
pos_esnps_adj_ea_z <- pos_esnps_results_data$adj_ea_z
names(pos_esnps_adj_ea_z) <- pos_esnps_results_data$ld_locus_id
pos_only_adj_ea_z <- pos_only_results_data$adj_ea_z
names(pos_only_adj_ea_z) <- pos_only_results_data$ld_locus_id

# Perform the adjusted GSEA for Positional + eSNPs
set.seed(1992)
pos_esnps_ea_adj_gsea <- fgsea(pos_esnps_target_sets,
                               pos_esnps_adj_ea_z, nperm = 10000,
                               scoreType = "pos")

# Now for Positional
set.seed(2004)
pos_only_ea_adj_gsea <- fgsea(pos_only_target_sets,
                              pos_only_adj_ea_z, nperm = 10000,
                              scoreType = "pos")

# Save tables of results:
pos_esnps_ea_adj_gsea %>%
  make_gsea_table %>%
  mutate(phenotype = "EA", positional_type = "Positional + eSNPs") %>%
  dplyr::select(phenotype, positional_type, everything()) %>%
  write_csv("data/gsea/ea_positional_esnps_results.csv")

pos_only_ea_adj_gsea %>%
  make_gsea_table %>%
  mutate(phenotype = "EA", positional_type = "Positional") %>%
  dplyr::select(phenotype, positional_type, everything()) %>%
  write_csv("data/gsea/ea_positional_results.csv")

# Now create a plot with the -log10 of the adjusted p-values:
adj_ea_gsea_results <- pos_esnps_ea_adj_gsea %>%
  simplify_gsea_results() %>%
  mutate(assign_type = "Positional + eSNPs") %>%
  bind_rows(
    pos_only_ea_adj_gsea %>%
      simplify_gsea_results() %>%
      mutate(assign_type = "Positional")
  ) %>%
  group_by(pathway) %>%
  mutate(total_logpval = sum(-log10(padj))) %>%
  ungroup() %>%
  mutate(pathway = fct_reorder(pathway, total_logpval)) %>%
  ggplot(aes(x = pathway, y = -log10(padj),
             color = assign_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = .75),
           width = 0.1, aes(fill = assign_type)) +
  geom_point(position = position_dodge(width = .75),
             size = 6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray", size = 1.5) +
  labs(y = "-log10(BH adjusted p-value)",
       x = "Gene/locus set") +
  scale_fill_manual(values = c("darkred", "darkblue"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = c("darkred", "darkblue"), guide = FALSE) +
  coord_flip() +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(.85, .15),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

# Create enrichment score plots -------------------------------------------

# Main display will include synaptic and epigenetic for Positional + eSNPs
# then supplement will have Positional plots

pos_esnps_synaptic_plot <- plotEnrichment(pos_esnps_target_sets[["synaptic"]],
                                          pos_esnps_adj_asd_z) +
  labs(x = "Adjusted ASD z statistic rank (in descending order)",
       y = "Synaptic enrichment score")

pos_esnps_epigenetic_plot <- plotEnrichment(pos_esnps_target_sets[["epigenetic"]],
                                            pos_esnps_adj_asd_z) +
  labs(x = "Adjusted ASD z statistic rank (in descending order)",
       y = "Epigenetic enrichment score")

# Combine with the enrichment pval plot for main display:

asd_gsea_plot <-
  plot_grid(adj_asd_gsea_results +
              guides(fill = guide_legend(
                reverse = TRUE,
                keywidth  = .5,
                keyheight = .5,
                default.unit = "inch")) +
              theme(legend.title = element_blank(),
                    legend.position = c(.7, .2),
                    legend.text = element_text(size = 18),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 16)),
            plot_grid(pos_esnps_synaptic_plot +
                        theme(axis.title.x = element_blank(),
                              axis.text.x = element_blank(),
                              axis.ticks.x = element_blank()),
                      NULL,
                      pos_esnps_epigenetic_plot ,
                      ncol = 1, rel_heights = c(1, -.1, 1),
                      align = "hv", labels = c("B", "", "C"),
                      label_fontface = "plain",
                      label_size = 18),
            ncol = 2, rel_widths = c(1.5, 1), label_size = 18,
            labels = c("A", ""), label_fontface = "plain")

save_plot("figures/main/f6_asd_gsea_results.jpg",
          asd_gsea_plot, ncol = 2, nrow = 2, base_asp = 2.5)
save_plot("figures/main/f6_asd_gsea_results.pdf",
          asd_gsea_plot, ncol = 2, nrow = 2, base_asp = 2.5)


# Create the positional versions for supplement:
pos_synaptic_plot <- plotEnrichment(pos_only_target_sets[["synaptic"]],
                                    pos_only_adj_asd_z) +
  labs(x = "Adjusted ASD z statistic rank (in descending order)",
       y = "Synaptic enrichment score")

pos_epigenetic_plot <- plotEnrichment(pos_only_target_sets[["epigenetic"]],
                                      pos_only_adj_asd_z) +
  labs(x = "Adjusted ASD z statistic rank (in descending order)",
       y = "Epigenetic enrichment score")

# Combine with the enrichment pval plot for main display:

asd_pos_gsea_plot <-
  plot_grid(pos_esnps_synaptic_plot +
              theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank()),
            NULL,
            pos_esnps_epigenetic_plot ,
            ncol = 1, rel_heights = c(1, -.1, 1),
            align = "hv", labels = c("A", "", "B"),
            label_fontface = "plain",
            label_size = 18)

save_plot("figures/suppl/f_pos_asd_gsea.jpg",
          asd_pos_gsea_plot, ncol = 1, nrow = 2, base_asp = 2.5)
save_plot("figures/suppl/f_pos_asd_gsea.pdf",
          asd_pos_gsea_plot, ncol = 1, nrow = 2, base_asp = 2.5)

