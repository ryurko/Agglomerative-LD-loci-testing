# PURPOSE: Generate the list of LD induced loci for each chromosome, based on
#          the Positional + eSNPs assignment

# SERVER EDITION

library(tidyverse)
library(data.table)
library(parallel)
library(snpcombineR)

# Source the helper functions
source("code/s4_create_ld_induced_genes_loci/agglom_algo_functions.R")

flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}


# Loop through chr to generate the separate results -----------------------

chr_loop <-
  mclapply(1:22, mc.cores = 22,
           function(chr_i) {

             # First load the tidy SNP-gene data for the chr
             tidy_chr_snp_gene_data <-
               read_csv(
                 paste0("data/tidy_snp_gene/positional_esnps/chr",
                        chr_i, "_snp_gene.csv"),
                 col_types = list(col_character(), col_character(),
                                  col_double(), col_double()))

             # Get the genotype data for the chromosome
             numeric_snp_ref_genotypes <- as.matrix(
               data.table::fread(
                 paste0("data/ld_reference/chr_level/chr",
                        chr_i, "_ref_genotypes.csv")))

             # Get the genotype data:
             numeric_snp_ref_genotypes <-
               numeric_snp_ref_genotypes[, unique(tidy_chr_snp_gene_data$snp_id)]
             numeric_snp_ref_genotypes <- flip_matrix(numeric_snp_ref_genotypes)
             
             # First get the chr gene correlation table (default window is 6 Mb)
             chr_gene_cor_table <- compute_inter_gene_cor(tidy_chr_snp_gene_data,
                                                          numeric_snp_ref_genotypes)
             
             # Save this file:
             write_csv(chr_gene_cor_table,
                       paste0(
                         "data/inter_gene_correlation/positional_esnps/chr",
                         chr_i, "_gene_cor.csv"))

             # Now construct the LD loci lists (using default thresholds and window)
             chr_ld_loci_results <-
               merge_correlated_genes(tidy_chr_snp_gene_data, numeric_snp_ref_genotypes,
                                      chr_gene_cor_table)

             saveRDS(chr_ld_loci_results,
                     paste0(
                       "data/merged_ld_loci/positional_esnps/chr",
                       chr_i, "_results.rds"))

             return(chr_i)

           })
