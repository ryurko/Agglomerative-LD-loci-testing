# PURPOSE: Initialize the tidy SNP-gene datasets using the two approaches -
#          (1) Positional and (2) Positional + eSNPs. Since the reference data
#          is based on hg19, each dataset will have the hg19 coordinates.

library(tidyverse)
library(data.table)

# Load the GWAS datasets --------------------------------------------------

# First the original set but remove SNPs with missing hg38 positions:
all_gwas_data <- fread("data/gwas/scz_asd_ea_gwas_results.csv") %>%
  .[hg38_id != "",]

# BrainVar eSNP-gene-data:
brainvar_esnp_gene_data <- fread("data/eqtls/brainvar/esnp_eqtl_results.csv")

# GTEx eSNP-gene pairs:
gtex_fcba9_esnp_gene_data <- fread("data/eqtls/gtex/fcba9_esnp_eqtl_results.csv")
gtex_accba24_esnp_gene_data <- fread("data/eqtls/gtex/accba24_esnp_eqtl_results.csv")


# Load the GENCODE table --------------------------------------------------

gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Create positional SNP-gene data -----------------------------------------

# First assign the SNPs from the GWAS data to GENCODE v21 genes based on their
# hg38 positions

# Create separate columns for the hg38 chr and bp:
snp_hg38_id_cols <- all_gwas_data[, .(hg38_id)] %>%
  as_tibble() %>%
  separate(hg38_id, c("hg38_chr", "hg38_bp"), remove = FALSE, sep = ":") %>%
  separate(hg38_bp, c("hg38_bp", "hg38_bp2"), remove = TRUE, sep = "-") %>%
  dplyr::select(-hg38_bp2) %>%
  mutate(hg38_chr = str_remove(hg38_chr, "chr"),
         # Catch the odd marked CHR:
         hg38_chr = as.numeric(sapply(hg38_chr,
                                      function(x) unlist(str_split(
                                        x, "_"
                                      ))[1])),
         hg38_bp = as.numeric(hg38_bp))

# Create new data.table objects with the necessary columns for finding overlaps:
# (see example here: https://www.biostars.org/p/371059/)

# Next initialize the SNP dt object for using foverlaps
snp_dt <- data.table(
  Chr = snp_hg38_id_cols$hg38_chr,
  start = snp_hg38_id_cols$hg38_bp,
  end = snp_hg38_id_cols$hg38_bp,
  id = snp_hg38_id_cols$hg38_id,
  key = c("Chr", "start", "end")
)

# Set-up the gene tables using the gencode data:
genes_dt <- data.table(
  Chr = gene_type_table$gene_chr,
  start = gene_type_table$start,
  end = gene_type_table$end,
  id = gene_type_table$ensembl_id,
  key = c("Chr", "start", "end")
)


# Make the overlaps table
snp_gene_assign <- foverlaps(snp_dt, genes_dt)

# Convert to the tidy SNP-gene dataset:
tidy_positional_snp_gene_data <- snp_gene_assign %>%
  # Remove SNPs unable to map to genes:
  .[!is.na(id),] %>%
  as_tibble() %>%
  dplyr::select(id, i.id, Chr) %>%
  dplyr::rename(ensembl_id = id,
                hg38_id = i.id,
                hg38_chr = Chr) %>%
  # Join the summary statistics:
  left_join(dplyr::select(all_gwas_data, hg38_id, chr,
                          hg19_id, scz_p, scz_or, scz_se,
                          asd_p, asd_or, asd_se, ea_p,
                          ea_beta, ea_se),
            by = "hg38_id")
length(unique(tidy_positional_snp_gene_data$hg38_id))
# [1] 2779780
length(unique(tidy_positional_snp_gene_data$hg19_id))
# [1] 2779780
length(unique(tidy_positional_snp_gene_data$ensembl_id))
# [1] 40581

# Save this:
write_csv(tidy_positional_snp_gene_data,
          "data/tidy_snp_gene/positional.csv")


# Next make tidy eSNP-gene datasets ---------------------------------------

# First for BrainVar:
tidy_brainvar_esnp_gene_data <- brainvar_esnp_gene_data %>%
  as_tibble() %>%
  dplyr::select(ensembl_id, esnp_hg38_id) %>%
  rename(hg38_id = esnp_hg38_id) %>%
  separate(hg38_id, c("hg38_chr", "hg38_bp"), remove = FALSE, sep = ":") %>%
  mutate(hg38_chr = str_remove(hg38_chr, "chr"),
         # Catch the odd marked CHR:
         hg38_chr = as.numeric(sapply(hg38_chr,
                                      function(x) unlist(str_split(
                                        x, "_"
                                      ))[1]))) %>%
  dplyr::select(-hg38_bp) %>%
  # Join the summary statistics:
  left_join(dplyr::select(all_gwas_data, hg38_id, chr,
                          hg19_id, scz_p, scz_or, scz_se,
                          asd_p, asd_or, asd_se, ea_p,
                          ea_beta, ea_se),
            by = "hg38_id")
length(unique(tidy_brainvar_esnp_gene_data$hg38_id))
# [1] 123664
length(unique(tidy_brainvar_esnp_gene_data$hg19_id))
# [1] 123664
length(unique(tidy_brainvar_esnp_gene_data$ensembl_id))
# [1] 6660

# Save this:
write_csv(tidy_brainvar_esnp_gene_data,
          "data/tidy_snp_gene/brainvar_eqtls.csv")

# Next for GTEx - make them separately for each tissue first before combining

# First fcba9:
tidy_gtex_fcba9_esnp_gene_data <- gtex_fcba9_esnp_gene_data %>%
  as_tibble() %>%
  dplyr::select(ensembl_id, hg19_id) %>%
  # Join the summary statistics:
  left_join(dplyr::select(all_gwas_data, hg38_id, chr,
                          hg19_id, scz_p, scz_or, scz_se,
                          asd_p, asd_or, asd_se, ea_p,
                          ea_beta, ea_se),
            by = "hg19_id") %>%
  separate(hg38_id, c("hg38_chr", "hg38_bp"), remove = FALSE, sep = ":") %>%
  mutate(hg38_chr = str_remove(hg38_chr, "chr"),
         # Catch the odd marked CHR:
         hg38_chr = as.numeric(sapply(hg38_chr,
                                      function(x) unlist(str_split(
                                        x, "_"
                                      ))[1]))) %>%
  dplyr::select(-hg38_bp)
# Save this:
write_csv(tidy_gtex_fcba9_esnp_gene_data,
          "data/tidy_snp_gene/gtex_fcba9_eqtls.csv")

# Next for accba24:
tidy_gtex_accba24_esnp_gene_data <- gtex_accba24_esnp_gene_data %>%
  as_tibble() %>%
  dplyr::select(ensembl_id, hg19_id) %>%
  # Join the summary statistics:
  left_join(dplyr::select(all_gwas_data, hg38_id, chr,
                          hg19_id, scz_p, scz_or, scz_se,
                          asd_p, asd_or, asd_se, ea_p,
                          ea_beta, ea_se),
            by = "hg19_id") %>%
  separate(hg38_id, c("hg38_chr", "hg38_bp"), remove = FALSE, sep = ":") %>%
  mutate(hg38_chr = str_remove(hg38_chr, "chr"),
         # Catch the odd marked CHR:
         hg38_chr = as.numeric(sapply(hg38_chr,
                                      function(x) unlist(str_split(
                                        x, "_"
                                      ))[1]))) %>%
  dplyr::select(-hg38_bp)
# Save this:
write_csv(tidy_gtex_accba24_esnp_gene_data,
          "data/tidy_snp_gene/gtex_accba24_eqtls.csv")

# Stack together for combined GTEx data:
tidy_gtex_esnp_gene_data <- tidy_gtex_fcba9_esnp_gene_data %>%
  bind_rows(tidy_gtex_accba24_esnp_gene_data) %>%
  distinct()

# Save:
write_csv(tidy_gtex_esnp_gene_data,
          "data/tidy_snp_gene/gtex_eqtls.csv")

# Combine the BrainVar and GTEx datasets:
tidy_esnp_gene_data <- tidy_brainvar_esnp_gene_data %>%
  bind_rows(tidy_gtex_esnp_gene_data) %>%
  distinct()

# Save:
write_csv(tidy_esnp_gene_data,
          "data/tidy_snp_gene/esnps.csv")


# Combine positional and eSNPs together -----------------------------------

tidy_positional_esnp_gene_data <- tidy_positional_snp_gene_data %>%
  bind_rows(tidy_esnp_gene_data) %>%
  distinct()

# And save:
write_csv(tidy_positional_esnp_gene_data,
          "data/tidy_snp_gene/positional_esnps.csv")





