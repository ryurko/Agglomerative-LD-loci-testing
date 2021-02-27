# PURPOSE: Initialize the SNP-gene input tables for generating the induced
#          LD gene-sets by modifying the column names of the current tables

library(tidyverse)


# Load the tidy SNP-gene data ---------------------------------------------

tidy_positional_snp_gene_data <-
  read_csv("data/tidy_snp_gene/positional.csv")

tidy_positional_esnp_gene_data <-
  read_csv("data/tidy_snp_gene/positional_esnps.csv")


# Load the GENCODE table --------------------------------------------------

gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Select and rename necessary columns -------------------------------------

positional_server_table <- tidy_positional_snp_gene_data %>%
  dplyr::select(ensembl_id, hg19_id, chr) %>%
  # Join the gene start and end:
  dplyr::left_join(dplyr::select(gene_type_table, ensembl_id, start, end),
                   by = "ensembl_id") %>%
  dplyr::rename(snp_id = hg19_id, gene_chr = chr,
                gene_id = ensembl_id)

positional_esnps_server_table <- tidy_positional_esnp_gene_data %>%
  dplyr::select(ensembl_id, hg19_id, chr) %>%
  # Join the gene start and end:
  dplyr::left_join(dplyr::select(gene_type_table, ensembl_id, start, end),
                   by = "ensembl_id") %>%
  dplyr::rename(snp_id = hg19_id, gene_chr = chr,
                gene_id = ensembl_id)

# Save both:
write_csv(positional_server_table,
          "data/tidy_snp_gene/server_version/positional.csv")
write_csv(positional_esnps_server_table,
          "data/tidy_snp_gene/server_version/positional_esnps.csv")

