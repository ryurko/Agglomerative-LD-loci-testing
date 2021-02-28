# PURPOSE: Initialize the LD reference data for the considered SNPs

library(tidyverse)
library(data.table)
library(snpStats)


# Load the tidy positional + eSNPs gene data ------------------------------
tidy_positional_esnp_gene_data <-
  read_csv("data/tidy_snp_gene/positional_esnps.csv")

# Initialize chromosome level genotype data -------------------------------

g1000_eur_fam_path <- "data/g1000_eur/g1000_eur.fam"
g1000_eur_bim_path <- "data/g1000_eur/g1000_eur.bim"
g1000_eur_bed_path <- "data/g1000_eur/g1000_eur.bed"

# Read in the PLINK data using snpStats:
ref_snps_plink <- read.plink(g1000_eur_bed_path,
                             g1000_eur_bim_path,
                             g1000_eur_fam_path)

# Obtain the SnpMatrix object (genotypes) table from ref_snps_plink list
ref_genotypes <- ref_snps_plink$genotypes

# Obtain the SNP information from ref_snps_plink list
ref_genobim <- ref_snps_plink$map
rm(ref_snps_plink)
colnames(ref_genobim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
ref_genobim <- as.data.table(ref_genobim)
ref_genobim[, hg19_id := paste0("chr", chr, ":", position, "-", position)]
all(unique(tidy_positional_esnp_gene_data$hg19_id) %in% ref_genobim$hg19_id)
# [1] TRUE

# Filter ref_genobim to the SNPs across the union of adult and fetal:
snp_ref_genobim <- ref_genobim[hg19_id %in%
                                 unique(tidy_positional_esnp_gene_data$hg19_id),]

all(snp_ref_genobim$SNP %in% colnames(ref_genotypes))
# [1] TRUE
rm(ref_genobim)

# Filter the genotype data to be these SNPs:
snp_ref_genotypes <- ref_genotypes[, snp_ref_genobim$SNP]
# Remove the large version:
rm(ref_genotypes)
dim(snp_ref_genotypes)
# [1]   503 2888822

# Create the numeric version of this matrix with alleles encoded as 0, 1, or 2:
numeric_snp_ref_genotypes <- as(snp_ref_genotypes, "numeric")
rm(snp_ref_genotypes)

# Remove the row names:
rownames(numeric_snp_ref_genotypes) <- NULL
# Change the colnames to be the hg19_id instead for consistency:
colnames(numeric_snp_ref_genotypes) <- snp_ref_genobim$hg19_id
rm(snp_ref_genobim)

# Next find and remove problematic SNPs - where all the values are the same:
problem_snps_i <- which(apply(numeric_snp_ref_genotypes, 2,
                              function(x) length(unique(x))) == 1)
# None!

# Save this dataset:
write_csv(as_tibble(numeric_snp_ref_genotypes),
          "data/ld_reference/ref_genotypes.csv")

all(tidy_positional_esnp_gene_data$hg19_id %in% colnames(numeric_snp_ref_genotypes))
# [1] TRUE

# Now make one for each chromosome:
for (chr_i in 1:22) {
  
  chr_snps <- tidy_positional_esnp_gene_data %>%
    filter(chr == chr_i) %>%
    pull(hg19_id) %>%
    unique()
  
  write_csv(as_tibble(
    numeric_snp_ref_genotypes[, chr_snps]),
    paste0("data/ld_reference/chr_level/chr",
           chr_i, "_ref_genotypes.csv"))
}










