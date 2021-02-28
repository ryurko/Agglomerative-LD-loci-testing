# PURPOSE: Initialize the GTEx eQTL datasets with the eSNPs and associated
#          genes as detected by applying BH to all candidate SNP-gene pairs,
#          analogous to the selection step by Werling et al for BrainVar.

library(tidyverse)
library(data.table)

# Load GWAS data ----------------------------------------------------------

scz_asd_ea_gwas_data <- fread("data/gwas/scz_asd_ea_gwas_results.csv")
# Note that fread converts NA to "" upon reading in
any(scz_asd_ea_gwas_data$hg38_id == "")
# [1] TRUE
length(which(scz_asd_ea_gwas_data$hg38_id == ""))
# [1] 2880

# Load GTEx data and select eSNP-gene pairs with BH -----------------------

# GTEx v7 based on hg19 and GENCODE v19: https://www.nature.com/articles/nature24277

# First load the Frontal Cortex BA9 dataset:
gtex_fcba9_all_pairs_data <-
  fread("gunzip -c data/eqtls/gtex/Brain_Frontal_Cortex_BA9.allpairs.txt.gz")

# Generate the BH q-values
gtex_fcba9_all_pairs_data[, bh_qval := p.adjust(pval_nominal,
                                                method = "BH")]
# Return the BH discoveries:
gtex_fcba9_bh_pairs_data <- gtex_fcba9_all_pairs_data %>%
  .[bh_qval <= 0.05,]
# Now free up memory by removing the original dataset:
rm(gtex_fcba9_all_pairs_data)

# Next the Anterior cingulate cortex BA24 dataset:
gtex_accba24_all_pairs_data <-
  fread("gunzip -c data/eqtls/gtex/Brain_Anterior_cingulate_cortex_BA24.allpairs.txt.gz")

# Generate the BH q-values
gtex_accba24_all_pairs_data[, bh_qval := p.adjust(pval_nominal,
                                                  method = "BH")]
# Return the BH discoveries:
gtex_accba24_bh_pairs_data <- gtex_accba24_all_pairs_data %>%
  .[bh_qval <= 0.05,]
# Now free up memory by removing the original dataset:
rm(gtex_accba24_all_pairs_data)


# Match the GWAS and GTEx data eSNPs --------------------------------------

# Proceed to first match on the hg19_id and alleles, then use the already lifted
# hg38_ids from the GWAS dataset to filter on which SNPs are eligible to use

gtex_fcba9_bh_pairs_data <- gtex_fcba9_bh_pairs_data %>%
  as_tibble() %>%
  separate(variant_id, c("chr", "bp", "a1", "a2", "build"), sep = "_")
gtex_accba24_bh_pairs_data <- gtex_accba24_bh_pairs_data %>%
  as_tibble() %>%
  separate(variant_id, c("chr", "bp", "a1", "a2", "build"), sep = "_")

# Create the hg19_id columns for these dataset:
gtex_fcba9_bh_pairs_data <- gtex_fcba9_bh_pairs_data %>%
  mutate(hg19_id = paste0("chr", chr, ":", bp, "-", bp))
gtex_accba24_bh_pairs_data <- gtex_accba24_bh_pairs_data %>%
  mutate(hg19_id = paste0("chr", chr, ":", bp, "-", bp))

# Add complement allele columns:
allele_complements <- c("A", "T", "C", "G")
names(allele_complements) <- c("T", "A", "G", "C")
gtex_fcba9_bh_pairs_data <- gtex_fcba9_bh_pairs_data %>%
  mutate(ca1 = as.character(allele_complements[a1]),
         ca2 = as.character(allele_complements[a2]))
gtex_accba24_bh_pairs_data <- gtex_accba24_bh_pairs_data %>%
  mutate(ca1 = as.character(allele_complements[a1]),
         ca2 = as.character(allele_complements[a2]))

# Grab the unique variant positions across these two datasets:
gtex_bh_variant_ids <- unique(c(gtex_fcba9_bh_pairs_data$hg19_id,
                                gtex_accba24_bh_pairs_data$hg19_id))
length(gtex_bh_variant_ids)
# [1] 474912
# Large overlap in eSNPs

# Next which GWAS SNPs are in the GTEx eSNPs based on the hg19 ids
gtex_esnp_gwas_data <- scz_asd_ea_gwas_data %>%
  .[hg19_id %in% gtex_bh_variant_ids,]
nrow(gtex_esnp_gwas_data)
# [1] 317286

any(gtex_esnp_gwas_data$hg38_id == "")
# [1] TRUE
length(which(gtex_esnp_gwas_data$hg38_id == ""))
# [1] 365

# Leave these SNPs in for now - to first match on the alleles prior to then
# removing any SNPs that do NOT liftover to hg38

# Now filter the gtex datasets to only be the variants in this data
# (note the only reason we are doing this filtering afterwards is to address
# the issue with being more comparable to BrainVar which provided us with
# the cis-eQTL BH pairs, rather than us making our own cutoff)
gtex_fcba9_bh_pairs_data <- gtex_fcba9_bh_pairs_data %>%
  filter(hg19_id %in% gtex_esnp_gwas_data$hg19_id)
gtex_accba24_bh_pairs_data <- gtex_accba24_bh_pairs_data %>%
  filter(hg19_id %in% gtex_esnp_gwas_data$hg19_id)

# Next - need to match on the alleles, can do this by joining the ASD alleles
# to each of the datasets:
gtex_fcba9_bh_pairs_data <- gtex_fcba9_bh_pairs_data %>%
  left_join(dplyr::select(gtex_esnp_gwas_data,
                          hg19_id, asd_a1, asd_a2) %>%
              rename(gwas_a1 = asd_a1,
                     gwas_a2 = asd_a2),
            by = "hg19_id") %>%
  # Create the match allele columns:
  mutate(match_allele_o = as.numeric(((a1 == gwas_a1) & (a2 == gwas_a2)) |
                                       ((a1 == gwas_a2) & (a2 == gwas_a1))),
         match_allele_c = as.numeric(((ca1 == gwas_a1) & (ca2 == gwas_a2)) |
                                       ((ca1 == gwas_a2) & (ca2 == gwas_a1))))
table(gtex_fcba9_bh_pairs_data$match_allele_o)
#   0      1
#   5 324298
table(gtex_fcba9_bh_pairs_data$match_allele_c)
#      0      1
# 276782  47521
with(gtex_fcba9_bh_pairs_data, table(match_allele_o, match_allele_c))
#               match_allele_c
# match_allele_o      0      1
#              0      5      0
#              1 276777  47521
# No impact! Can just use match_allele_o

gtex_fcba9_match_gwas_data <- gtex_fcba9_bh_pairs_data %>%
  filter(match_allele_o == 1)

# Repeat for the other set:
gtex_accba24_bh_pairs_data <- gtex_accba24_bh_pairs_data %>%
  left_join(dplyr::select(gtex_esnp_gwas_data,
                          hg19_id, asd_a1, asd_a2) %>%
              rename(gwas_a1 = asd_a1,
                     gwas_a2 = asd_a2),
            by = "hg19_id") %>%
  # Create the match allele columns:
  mutate(match_allele_o = as.numeric(((a1 == gwas_a1) & (a2 == gwas_a2)) |
                                       ((a1 == gwas_a2) & (a2 == gwas_a1))),
         match_allele_c = as.numeric(((ca1 == gwas_a1) & (ca2 == gwas_a2)) |
                                       ((ca1 == gwas_a2) & (ca2 == gwas_a1))))
table(gtex_accba24_bh_pairs_data$match_allele_o)
#   0      1
#   2 265468
table(gtex_accba24_bh_pairs_data$match_allele_c)
#      0      1
# 226449  39021
with(gtex_accba24_bh_pairs_data, table(match_allele_o, match_allele_c))
#               match_allele_c
# match_allele_o      0      1
#              0      2      0
#              1 226447  39021
# No impact! Can just use match_allele_o
gtex_accba24_match_gwas_data <- gtex_accba24_bh_pairs_data %>%
  filter(match_allele_o == 1)

# Now filter the GWAS data only to the SNPs with matching alleles:
match_gtex_bh_variant_ids <- unique(c(gtex_fcba9_match_gwas_data$hg19_id,
                                      gtex_accba24_match_gwas_data$hg19_id))
length(match_gtex_bh_variant_ids)
# [1] 317280

gtex_esnp_gwas_data <- gtex_esnp_gwas_data %>%
  .[hg19_id %in% match_gtex_bh_variant_ids,]
nrow(gtex_esnp_gwas_data)
# [1] 317280

any(gtex_esnp_gwas_data$hg38_id == "")
# [1] TRUE
length(which(gtex_esnp_gwas_data$hg38_id == ""))
# [1] 365
# Same count as before post allele matching

# Remove these SNPs
gtex_esnp_gwas_data <- gtex_esnp_gwas_data %>%
  filter(hg38_id != "")
nrow(gtex_esnp_gwas_data)
# [1] 316915

# Update the eQTL gene pair datasets accordingly as well:
gtex_fcba9_match_gwas_data <- gtex_fcba9_match_gwas_data %>%
  filter(hg19_id %in% gtex_esnp_gwas_data$hg19_id)
gtex_accba24_match_gwas_data <- gtex_accba24_match_gwas_data %>%
  filter(hg19_id %in% gtex_esnp_gwas_data$hg19_id)


# Check eQTLs with GENCODE genes ------------------------------------------

# Load the GENCODE gene table:
gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Next need to convert GTEx gene ids to the base Ensembl id:
gtex_fcba9_match_gwas_data <- gtex_fcba9_match_gwas_data %>%
  mutate(ensembl_id = str_sub(gene_id, 1, 15))
gtex_accba24_match_gwas_data <- gtex_accba24_match_gwas_data %>%
  mutate(ensembl_id = str_sub(gene_id, 1, 15))

# Check to see which of the eQTL genes are not in the GENCODE v21 table -
# and for ease will just drop those from consideration:
all(unique(gtex_fcba9_match_gwas_data$ensembl_id) %in% gene_type_table$ensembl_id)
# [1] FALSE
all(unique(gtex_accba24_match_gwas_data$ensembl_id) %in% gene_type_table$ensembl_id)
# [1] FALSE

# Both false - so now need to remove the missing genes:
gtex_eqtl_genes <- unique(c(gtex_fcba9_match_gwas_data$ensembl_id,
                            gtex_accba24_match_gwas_data$ensembl_id))
length(gtex_eqtl_genes)
# [1] 9199

# Missing from GENCODE v21:
missing_gtex_eqtl_genes <- setdiff(gtex_eqtl_genes, gene_type_table$ensembl_id)
length(missing_gtex_eqtl_genes)
# [1] 187

# Next check biomaRt to see if any of these genes have hg38 coordinates:
library(biomaRt)
ensembl_mart_hg38 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_genes_hg38 <- getBM(attributes = c("ensembl_gene_id",
                                           "chromosome_name",
                                           "start_position", "end_position"),
                            filters = "ensembl_gene_id",
                            values = missing_gtex_eqtl_genes,
                            mart = ensembl_mart_hg38)
nrow(ensembl_genes_hg38)
# [1] 0 - nothing! so these genes can be removed

# Will drop these from the datasets:
clean_gtex_fcba9_gwas_data <- gtex_fcba9_match_gwas_data %>%
  filter(ensembl_id %in% gene_type_table$ensembl_id)
clean_gtex_accba24_gwas_data <- gtex_accba24_match_gwas_data %>%
  filter(ensembl_id %in% gene_type_table$ensembl_id)

# Now only consider the remaining SNPs:
clean_gtex_bh_variant_ids <- unique(c(clean_gtex_fcba9_gwas_data$hg19_id,
                                      clean_gtex_accba24_gwas_data$hg19_id))
length(clean_gtex_bh_variant_ids)
# [1] 313316, change of 317280 - 313316 = 3964 SNPs

clean_gtex_esnp_gwas_data <- gtex_esnp_gwas_data %>%
  filter(hg19_id %in% clean_gtex_bh_variant_ids)
nrow(clean_gtex_esnp_gwas_data)
# [1] 313316

# Save these datasets -----------------------------------------------------

write_csv(clean_gtex_esnp_gwas_data,
          "data/eqtls/gtex/esnp_gwas_data.csv")
write_csv(clean_gtex_fcba9_gwas_data,
          "data/eqtls/gtex/fcba9_esnp_eqtl_results.csv")
write_csv(clean_gtex_accba24_gwas_data,
          "data/eqtls/gtex/accba24_esnp_eqtl_results.csv")








