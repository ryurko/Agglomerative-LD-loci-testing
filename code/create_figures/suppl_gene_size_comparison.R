# PURPOSE: Compare the gene size distributions between selection and null results.
#          Then initialize the matching sets for FUMA input to adjust for gene size.

library(tidyverse)

# Load the tidy SNP gene and GENCODE tables -------------------------------

# First the tidy SNP-to-gene tables which will be used for computing gene-size
# in terms of the number of assigned SNPs instead of just length
tidy_snp_gene_data <-
  read_csv("data/tidy_snp_gene/positional_esnps.csv")

# Next the gene to LD loci data, to exclude any genes in selected LD loci:
tidy_gene_ld_loci_data <-
  read_csv("data/tidy_gene_ld_loci/positional_esnps/agglom_rsquared25_ld_loci.csv")

# Load the GENCODE data
gencode_data <- read_csv("data/gencode/gencode_v21_table.csv")

# Load the results for the selected LD loci -----------------------------

# ASD
asd_gene_ld_loci_results <-
  read_csv("data/adapt_results/positional_esnps/snp_gene_tables/asd_rsquared25_gene_ld_loci_table.csv")

# SCZ
scz_gene_ld_loci_results <-
  read_csv("data/adapt_results/positional_esnps/snp_gene_tables/scz_rsquared25_gene_ld_loci_table.csv")


# Load signal genes -------------------------------------------------------

asd_signals_genes <-
  read_lines("data/adapt_results/gene_lists/signals/positional_esnps/asd_rsquared25_disc_genes.txt")
scz_signals_genes <-
  read_lines("data/adapt_results/gene_lists/signals/positional_esnps/scz_rsquared25_disc_genes.txt")


# Construct tables of sizes by length and SNPs  ---------------------------

# First count the number of SNPs assigned to each gene:
gene_info_data <- tidy_snp_gene_data %>%
  group_by(ensembl_id) %>%
  summarize(n_snps = n()) %>%
  ungroup() %>%
  # join the LD loci id
  dplyr::left_join(tidy_gene_ld_loci_data, by = "ensembl_id") %>%
  # join the gencode info
  dplyr::left_join(gencode_data, by = "ensembl_id") %>%
  # Compute size based on end - start
  mutate(gene_size = end - start,
         # Make an indicator denoting the LD loci selection
         is_asd_ld_loci = ifelse(ld_loci_id %in% asd_gene_ld_loci_results$ld_loci_id,
                            "Yes", "No"),
         is_scz_ld_loci = ifelse(ld_loci_id %in% scz_gene_ld_loci_results$ld_loci_id,
                            "Yes", "No"),
         # Indicators for signal genes
         is_asd_sig = ifelse(ensembl_id %in% asd_signals_genes, "Yes", "No"),
         is_scz_sig = ifelse(ensembl_id %in% scz_signals_genes, "Yes", "No"))

# Convert Ensembl to Entrez ids -------------------------------------------

unique_entrez_ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                           keys = gene_info_data$ensembl_id,
                                           column = "ENTREZID",
                                           keytype = "ENSEMBL",
                                           # Only use ids with a unique identifier
                                           multiVals = "asNA")
gene_info_data$entrez_id <- unique_entrez_ids



# Load GENCODE data -------------------------------------------------------

gencode_data <-
  read_csv("data/gencode/gencode_v21_table.csv")

# How many of these missing do not have unique Entrez gene IDs:
gencode_data$entrez_id <-
  AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                        keys = gencode_data$ensembl_id,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        # Only use ids with a unique identifier
                        multiVals = "asNA")

# Make the figure comparing the selected and null gene sizes --------------

# Compare the ECDFs between the genes in the selected LD loci and the
# "nulls" by the number of SNPs for both phenotypes. Need to stack the two
# long versions of the dataset together with the separate phenotype results
gene_info_data %>%
  dplyr::select(ensembl_id, n_snps, is_asd_ld_loci, is_scz_ld_loci) %>%
  pivot_longer(is_asd_ld_loci:is_scz_ld_loci,
               names_to = "phenotype", values_to = "is_selected") %>%
  mutate(phenotype = toupper(str_sub(phenotype, start = 4, end = 6))) %>%
  ggplot(aes(x = n_snps, color = is_selected)) +
  stat_ecdf() +
  facet_wrap(~ phenotype, ncol = 2, scales = "free_x") +
  theme_bw() +
  ggthemes::scale_color_colorblind() +
  labs(x = "Gene size", y = "Proportion of genes",
       color = "In selected LD loci?") +
  theme(legend.position = "bottom",
        strip.background = element_blank())

# Next with the signals
size_sig_ecdf <- bind_rows(
  gene_info_data %>%
    dplyr::select(ensembl_id, n_snps, is_asd_ld_loci, is_asd_sig) %>%
    filter(is_asd_ld_loci == "No" | is_asd_sig == "Yes") %>%
    mutate(gene_type = ifelse(is_asd_ld_loci == "No", "Non-implicated", "Signal genes"),
           phenotype = "ASD") %>%
    dplyr::select(-is_asd_ld_loci, -is_asd_sig),
  gene_info_data %>%
    dplyr::select(ensembl_id, n_snps, is_scz_ld_loci, is_scz_sig) %>%
    filter(is_scz_ld_loci == "No" | is_scz_sig == "Yes") %>%
    mutate(gene_type = ifelse(is_scz_ld_loci == "No", "Non-implicated", "Signal genes"),
           phenotype = "SCZ") %>%
    dplyr::select(-is_scz_ld_loci, -is_scz_sig)
) %>%
  ggplot(aes(x = n_snps, color = gene_type)) +
  stat_ecdf() +
  facet_wrap(~ phenotype, ncol = 2, scales = "free_x") +
  theme_bw() +
  ggthemes::scale_color_colorblind() +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(x = "Gene size", y = "Proportion of genes") +
  theme(legend.position = c(.2, .4),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_blank())
cowplot::save_plot("figures/suppl/f_gene_size_comparison.jpg",
                   size_sig_ecdf, ncol = 2, base_height = 3)
cowplot::save_plot("figures/suppl/f_gene_size_comparison.pdf",
                   size_sig_ecdf, ncol = 2, base_height = 3)
