# PURPOSE: Create and save the LD induced loci using Positional assignment

# WARNING: This code was used on a server distributing the construction of LD
#          induced loci for each chr in parallel. 

library(tidyverse)
library(data.table)
library(parallel)
library(snpcombineR)

# Source the main functions
source("code/s05_create_ld_induced_loci/form_ld_induced_loci.R")

# Init helper function for 
flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Load the dataset
tidy_snp_data <-
  read_csv("data/input/tidy_snp_gene/server_version/positional.csv",
           col_types = list(col_character(), col_character(), col_double(),
                            col_double(), col_double()))

# Will loop through each chromosome, computing and saving the gene correlations
# for all of the genes in each:
all_chr <- unique(tidy_snp_data$gene_chr)

# r^2 = 0.25 --------------------------------------------------------------

ld_loci_walk <-
  mclapply(1:22, mc.cores = 22,
           function(chr_i) {
             
             # Get the genotype data for the chromosome
             numeric_snp_ref_genotypes <- as.matrix(
               data.table::fread(
                 paste0("data/ld_reference/chr_level/chr", chr_i, "_ref_genotypes.csv")))
             
             
             # Get only this chr data:
             tidy_chr_snp_data <- tidy_snp_data %>%
               filter(gene_chr == chr_i) %>%
               dplyr::select(-gene_chr)
             
             # Get the genotype data:
             chr_genotype_data <-
               numeric_snp_ref_genotypes[, unique(tidy_chr_snp_data$snp_id)]
             chr_genotype_data <- flip_matrix(chr_genotype_data)
             
             # Form the LD induced loci
             chr_ld_loci_results <- 
               form_ld_induced_loci(tidy_chr_snp_data, chr_genotype_data,
                                    window_size = 6,
                                    # rho threshold of sqrt(.25)
                                    rho_threshold = sqrt(.25))
             # Save
             saveRDS(chr_ld_loci_results,
                     paste0("data/ld_induced_loci/rsquared25/positional/chr",
                            chr_i, "_results.rds"))
             
             return(chr_i)
             
           })


# r^2 = 0.50 --------------------------------------------------------------

ld_loci_walk <-
  mclapply(1:22, mc.cores = 22,
           function(chr_i) {
             
             # Get the genotype data for the chromosome
             numeric_snp_ref_genotypes <- as.matrix(
               data.table::fread(
                 paste0("data/ld_reference/chr_level/chr", chr_i, "_ref_genotypes.csv")))
             
             # Get only this chr data:
             tidy_chr_snp_data <- tidy_snp_data %>%
               filter(gene_chr == chr_i) %>%
               dplyr::select(-gene_chr)
             
             # Get the genotype data:
             chr_genotype_data <-
               numeric_snp_ref_genotypes[, unique(tidy_chr_snp_data$snp_id)]
             chr_genotype_data <- flip_matrix(chr_genotype_data)
             
             # Form the LD induced loci
             chr_ld_loci_results <- 
               form_ld_induced_loci(tidy_chr_snp_data, chr_genotype_data,
                                    window_size = 6,
                                    # rho threshold of sqrt(.5)
                                    rho_threshold = sqrt(.5))
             # Save
             saveRDS(chr_ld_loci_results,
                     paste0("data/ld_induced_loci/rsquared50/positional/chr",
                            chr_i, "_results.rds"))
             
             return(chr_i)
             
           })

# r^2 = 0.75 --------------------------------------------------------------

ld_loci_walk <-
  mclapply(1:22, mc.cores = 22,
           function(chr_i) {
             
             # Get the genotype data for the chromosome
             numeric_snp_ref_genotypes <- as.matrix(
               data.table::fread(
                 paste0("data/ld_reference/chr_level/chr", chr_i, "_ref_genotypes.csv")))
             
             # Get only this chr data:
             tidy_chr_snp_data <- tidy_snp_data %>%
               filter(gene_chr == chr_i) %>%
               dplyr::select(-gene_chr)
             
             # Get the genotype data:
             chr_genotype_data <-
               numeric_snp_ref_genotypes[, unique(tidy_chr_snp_data$snp_id)]
             chr_genotype_data <- flip_matrix(chr_genotype_data)
             
             # Form the LD induced loci
             chr_ld_loci_results <- 
               form_ld_induced_loci(tidy_chr_snp_data, chr_genotype_data,
                                    window_size = 6,
                                    # rho threshold of sqrt(.75)
                                    rho_threshold = sqrt(.75))
             # Save
             saveRDS(chr_ld_loci_results,
                     paste0("data/ld_induced_loci/rsquared75/positional/chr",
                            chr_i, "_results.rds"))
             
             return(chr_i)
             
           })

