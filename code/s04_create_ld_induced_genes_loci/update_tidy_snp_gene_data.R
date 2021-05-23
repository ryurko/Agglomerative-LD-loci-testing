# PURPOSE: Initialize the tidy SNP-gene datasets to use for merging on the server

library(tidyverse)

# Load the two versions based on positional assignment --------------------

tidy_positional_esnp_gene_data <-
  read_csv("data/tidy_snp_gene/positional_esnps.csv")

tidy_positional_gene_data <-
  read_csv("data/tidy_snp_gene/positional.csv")


# Load gencode data to join -----------------------------------------------

gencode_data <- read_csv("data/gencode/gencode_v21_table.csv")

# Create and save each tidy SNP-gene dataset separately -------------------

# Positional + eSNPs
tidy_snp_gene_positional_esnps <- tidy_positional_esnp_gene_data %>%
  dplyr::select(hg19_id, ensembl_id, chr) %>%
  # Join the gene start and end:
  dplyr::left_join(dplyr::select(gencode_data, ensembl_id, start, end),
                   by = "ensembl_id") %>%
  dplyr::rename(snp_id = hg19_id, gene_id = ensembl_id, gene_chr = chr)
write_csv(tidy_snp_gene_positional_esnps,
          "data/tidy_snp_gene/positional_esnps.csv")

# Repeat for Positional
tidy_snp_gene_positional <- tidy_positional_gene_data %>%
  dplyr::select(hg19_id, ensembl_id, chr) %>%
  # Join the gene start and end:
  dplyr::left_join(dplyr::select(gencode_data, ensembl_id, start, end),
                   by = "ensembl_id") %>%
  dplyr::rename(snp_id = hg19_id, gene_id = ensembl_id, gene_chr = chr)
write_csv(tidy_snp_gene_positional,
          "data/tidy_snp_gene/positional.csv")


# Save at the chromosome level --------------------------------------------

walk(1:22,
     function(chr_i) {

       # Drop the gene_chr column after filtering
       tidy_snp_gene_positional_esnps %>%
         filter(gene_chr == chr_i) %>%
         dplyr::select(-gene_chr) %>%
         write_csv(paste0("data/tidy_snp_gene/positional_esnps/chr",
                          chr_i, "_snp_gene.csv"))

       tidy_snp_gene_positional %>%
         filter(gene_chr == chr_i) %>%
         dplyr::select(-gene_chr) %>%
         write_csv(paste0("data/tidy_snp_gene/positional/chr",
                          chr_i, "_snp_gene.csv"))

     })

