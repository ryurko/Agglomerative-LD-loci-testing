# PURPOSE: Initialize the BrainVar eQTL datasets with the eSNPs and associated
#          genes reported by Werling et al in their public results.

library(tidyverse)
library(data.table)

# Load GWAS data ----------------------------------------------------------

scz_asd_ea_gwas_data <- fread("data/gwas/scz_asd_ea_gwas_results.csv")

# Load the public BrainVar data -------------------------------------------

# The following function is from code here:
# https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(sheet_i) readxl::read_excel(filename,
                                                           sheet = sheet_i))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}

# Load the different BrainVar sheets, first the eQTL sheets:
brainvar_eqtl_sheets <- read_excel_allsheets("data/eqtls/brainvar/eqtls.xlsx")

# Access separately the HCP eQTL sheet - containing all of the BrainVar eQTLs
# (meaning it is an eQTL regardless of the sample period)
brainvar_eqtl_data <- brainvar_eqtl_sheets$HCP_allEQTLs_FDR0.05

# Separate the GeneId into the ensembl_id and symbol columns, and also separte
# the Variant into esnp_chr, esnp_bp, esnp_a1, esnp_a2, then create the esnp_hg38_id
# and tss_hg38 columns to use for liftOver:
brainvar_eqtl_data <- brainvar_eqtl_data %>%
  separate(GeneId, c("ensembl_id", "symbol"), sep = "\\|") %>%
  separate(Variant, c("esnp_chr", "esnp_bp", "esnp_a1", "esnp_a2"), sep = ":",
           remove = FALSE) %>%
  mutate(esnp_chr = as.numeric(str_remove(esnp_chr, "chr")),
         esnp_hg38_id = paste0("chr", esnp_chr, ":", esnp_bp, "-", esnp_bp),
         tss_hg38 = paste0("chr", esnp_chr, ":", TSS, "-", TSS)) %>%
  # Remove missing chromosome information since those are sex chromosomes:
  filter(!is.na(esnp_chr))

# Add complement allele columns:
allele_complements <- c("A", "T", "C", "G")
names(allele_complements) <- c("T", "A", "G", "C")
brainvar_eqtl_data <- brainvar_eqtl_data %>%
  mutate(esnp_ca1 = as.character(allele_complements[esnp_a1]),
         esnp_ca2 = as.character(allele_complements[esnp_a2]))


# Match the GWAS and BrainVar data eSNPs ----------------------------------

# Next which GWAS SNPs are in the BrainVar eSNPs based on the hg38 ids
brainvar_esnp_gwas_data <- scz_asd_ea_gwas_data %>%
  .[hg38_id %in% brainvar_eqtl_data$esnp_hg38_id,]
nrow(brainvar_esnp_gwas_data)
# [1] 123666

# Now filter the brainvar to only be those that match these:
brainvar_eqtl_match_gwas_data <- brainvar_eqtl_data %>%
  filter(esnp_hg38_id %in% brainvar_esnp_gwas_data$hg38_id)

# Next go through the distinct matching BrainVar hg38_ids and alleles:
brainvar_eqtl_match_alleles <- brainvar_eqtl_match_gwas_data %>%
  dplyr::select(esnp_hg38_id, esnp_a1, esnp_a2, esnp_ca1, esnp_ca2) %>%
  distinct() %>%
  # join one set of the alleles (since they already match across asd, scz, ea)
  left_join(dplyr::select(brainvar_esnp_gwas_data, hg38_id, asd_a1, asd_a2),
            by = c("esnp_hg38_id" = "hg38_id")) %>%
  # Check allele matching in both ways:
  mutate(match_gwas_alleles_o = as.numeric(((esnp_a1 == asd_a1) & (esnp_a2 == asd_a2)) |
                                             ((esnp_a2 == asd_a1) & (esnp_a1 == asd_a2))),
         match_gwas_alleles_c = as.numeric(((esnp_ca1 == asd_a1) & (esnp_ca2 == asd_a2)) |
                                             ((esnp_ca2 == asd_a1) & (esnp_ca1 == asd_a2))))
table(brainvar_eqtl_match_alleles$match_gwas_alleles_o)
#   0      1
# 246 123420
table(brainvar_eqtl_match_alleles$match_gwas_alleles_c)
#      0      1
# 105250  18416
with(brainvar_eqtl_match_alleles, table(match_gwas_alleles_o, match_gwas_alleles_c))
#                      match_gwas_alleles_c
# match_gwas_alleles_o      0      1
#                    0      2    244
#                    1 105248  18172
# Ah so that so means there are some captured by this then!

# Filter to only the SNPs with the matching alleles:
brainvar_eqtl_match_alleles <- brainvar_eqtl_match_alleles %>%
  mutate(union_match_gwas_alleles = pmax(match_gwas_alleles_o, match_gwas_alleles_c)) %>%
  filter(union_match_gwas_alleles == 1)

# Now filter down to only these SNPs in the eSNP GWAS data:
brainvar_esnp_gwas_data <- brainvar_esnp_gwas_data %>%
  filter(hg38_id %in% brainvar_eqtl_match_alleles$esnp_hg38_id)

# And then in the eQTL data:
brainvar_eqtl_match_gwas_data <- brainvar_eqtl_match_gwas_data %>%
  filter(esnp_hg38_id %in% brainvar_eqtl_match_alleles$esnp_hg38_id)


# Check eQTLs with GENCODE genes ------------------------------------------

# Load the GENCODE gene table:
gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Are all of the BrainVar eQTL genes in this table?
all(unique(brainvar_eqtl_match_gwas_data$ensembl_id) %in%
      gene_type_table$ensembl_id)
# [1] TRUE
# All of the genes are in this table so no additional filtering is necessary

# Save these datasets -----------------------------------------------------

write_csv(brainvar_esnp_gwas_data,
          "data/eqtls/brainvar/esnp_gwas_data.csv")

write_csv(brainvar_eqtl_match_gwas_data,
          "data/eqtls/brainvar/esnp_eqtl_results.csv")





