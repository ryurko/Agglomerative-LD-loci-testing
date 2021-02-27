# PURPOSE: Initialize the dataset of GWAS results for ASD, SCZ, and EA, creating
#          a single dataset of SNPs with all three phenotypes. Filter SNPs to
#          only include those with MAF > 0.05 based on 1000 Genomes reference
#          data, since this will be the source of the LD structure.

# AUTHOR: Ron Yurko

library(tidyverse)
library(data.table)
library(snpStats)

# Load GWAS summary statistics --------------------------------------------

# ASD (hg19 - in readme)
asd_gwas_data <- fread("gunzip -c data/gwas/asd/iPSYCH-PGC_ASD_Nov2017.gz")

# SCZ (hg19 - in manuscript methods)
scz_gwas_data <- fread("gunzip -c data/gwas/scz/sczvscont-sumstat.gz")

# EA (hg19 - based on supplement referring to use of GENCODE v26 backmapped to GRCh37)
ea_gwas_data <- fread("data/gwas/ea/GWAS_EA_excl23andME.txt")


# Create SNP ids and merge datasets ---------------------------------------

# Create the hg19_id (chrCHR:BP-BP) for each:
asd_gwas_data[, hg19_id := paste0("chr", CHR, ":", BP, "-", BP)]
scz_gwas_data[, hg19_id := paste0("chr", CHR, ":", BP, "-", BP)]
ea_gwas_data[, hg19_id := paste0("chr", CHR, ":", POS, "-", POS)]

# What's the intersection of these SNPs:
overlap_hg19_ids <- intersect(asd_gwas_data$hg19_id,
                              intersect(scz_gwas_data$hg19_id,
                                        ea_gwas_data$hg19_id))

# Filter the datasets to only be these SNPs, and drop unnecessary columns:
asd_gwas_data <- asd_gwas_data[hg19_id %in% overlap_hg19_ids,
                               list(hg19_id, SNP, CHR, BP, A1, A2,
                                    OR, SE, P, INFO)]
scz_gwas_data <- scz_gwas_data[hg19_id %in% overlap_hg19_ids,
                               list(hg19_id, SNP, CHR, BP, A1, A2, OR, SE,
                                    P, INFO, FRQ_A_33426, FRQ_U_32541)]
ea_gwas_data <- ea_gwas_data[hg19_id %in% overlap_hg19_ids,
                             list(hg19_id, MarkerName, CHR, POS, A1, A2, Beta, SE,
                                  Pval, EAF)]
rm(overlap_hg19_ids)

# Rename columns:
setnames(asd_gwas_data,
         old = c('SNP', 'CHR', 'BP', 'A1', 'A2', 'OR', 'SE', 'P', 'INFO'),
         new = c('asd_rsid', 'chr', 'bp', 'asd_a1', 'asd_a2', 'asd_or',
                 'asd_se', 'asd_p', 'asd_info'))
setnames(scz_gwas_data,
         old = c('SNP', 'CHR', 'BP', 'A1', 'A2', 'OR', 'SE', 'P', 'INFO',
                 'FRQ_A_33426', 'FRQ_U_32541'),
         new = c('scz_rsid', 'chr', 'bp', 'scz_a1', 'scz_a2', 'scz_or',
                 'scz_se', 'scz_p', 'scz_info', 'scz_a1_case_freq',
                 'scz_a2_control_freq'))
setnames(ea_gwas_data,
         old = c('MarkerName', 'CHR', 'POS', 'A1', 'A2', 'Beta', 'SE', 'Pval', 'EAF'),
         new = c('ea_rsid', 'chr', 'bp', 'ea_a1', 'ea_a2', 'ea_beta',
                 'ea_se', 'ea_p', 'ea_a1_ref_freq'))

# Add complement allele columns:
allele_complements <- c("A", "T", "C", "G")
names(allele_complements) <- c("T", "A", "G", "C")

asd_gwas_data[,
              c("asd_ca1", "asd_ca2") :=
                .(as.character(allele_complements[asd_a1]),
                  as.character(allele_complements[asd_a2]))]
scz_gwas_data[,
              c("scz_ca1", "scz_ca2") :=
                .(as.character(allele_complements[scz_a1]),
                  as.character(allele_complements[scz_a2]))]
ea_gwas_data[,
             c("ea_ca1", "ea_ca2") :=
               .(as.character(allele_complements[ea_a1]),
                 as.character(allele_complements[ea_a2]))]

# Use these ids as keys
snp_keys <- c("hg19_id", "chr", "bp")
setkeyv(asd_gwas_data, snp_keys)
setkeyv(scz_gwas_data, snp_keys)
setkeyv(ea_gwas_data, snp_keys)

# Join the datasets together - useing left join starting with SCZ and ASD:
scz_asd_gwas_data <- merge(scz_gwas_data, asd_gwas_data,
                           all.x = TRUE)
# Filter to only be the SNPs with matching alleles in two ways:
# (1) either direction with given alleles, and (2) either direction with complementary
# (well technically its four ways)

# First observe how many many match in those two ways:
scz_asd_gwas_data[,
                  # First approach is the given alleles
                  match_scz_asd :=
                    as.numeric(((scz_a1 == asd_a1) &
                                  (scz_a2 == asd_a2)) |
                                 ((scz_a1 == asd_a2) &
                                    (scz_a2 == asd_a1)))]
table(scz_asd_gwas_data$match_scz_asd)
#    0       1
# 3102 6523119
# Next approach is the complementary alleles - for ASD:
scz_asd_gwas_data[,
                  # First approach is the given alleles
                  match_scz_asd_c :=
                    as.numeric(((scz_a1 == asd_ca1) &
                                  (scz_a2 == asd_ca2)) |
                                 ((scz_a1 == asd_ca2) &
                                    (scz_a2 == asd_ca1)))]
table(scz_asd_gwas_data$match_scz_asd_c)
#       0       1
# 5544674  981407
# Okay - so whats the union of the two?
scz_asd_gwas_data[,
                  union_match_scz_asd := pmax(match_scz_asd, match_scz_asd_c)]
table(scz_asd_gwas_data$union_match_scz_asd)
#    0       1
# 2962 6523119

# Upon inspection - the only missings are more than one base pair change (indels)
# so these can be safely ignored

# Filter to these matching alleles:
scz_asd_gwas_data <- scz_asd_gwas_data[union_match_scz_asd == 1,]
# Remove the asd and scz gwas data:
rm(scz_gwas_data, asd_gwas_data)

# Next join the EA SNPs and again filter on matching alleles for both scz and asd:
scz_asd_ea_gwas_data <- merge(scz_asd_gwas_data, ea_gwas_data,
                              all.x = TRUE) %>%
  # Do it both ways with both phenotypes:
  .[, ':='(match_scz_ea = as.numeric(((scz_a1 == ea_a1) &
                                        (scz_a2 == ea_a2)) |
                                       ((scz_a1 == ea_a2) &
                                          (scz_a2 == ea_a1))),
           match_scz_ea_c = as.numeric(((scz_a1 == ea_ca1) &
                                          (scz_a2 == ea_ca2)) |
                                         ((scz_a1 == ea_ca2) &
                                            (scz_a2 == ea_ca1))),
           match_asd_ea = as.numeric(((ea_a1 == asd_a1) &
                                        (ea_a2 == asd_a2)) |
                                       ((ea_a1 == asd_a2) &
                                          (ea_a2 == asd_a1))),
           match_asd_ea_c = as.numeric(((ea_ca1 == asd_a1) &
                                          (ea_ca2 == asd_a2)) |
                                         ((ea_ca1 == asd_a2) &
                                            (ea_ca2 == asd_a1))))]
table(scz_asd_ea_gwas_data$match_scz_ea)
#       1
# 6523119
table(scz_asd_ea_gwas_data$match_scz_ea_c)
# 0       1
# 5541712  981407
table(scz_asd_ea_gwas_data$match_asd_ea)
#       1
# 6523119
table(scz_asd_ea_gwas_data$match_asd_ea_c)
# 0       1
# 5541712  981407

# Okay so no difference at all!
scz_asd_ea_gwas_data[,
                     ':='(
                       union_match_scz_ea = pmax(match_scz_ea, match_scz_ea_c),
                       union_match_asd_ea = pmax(match_asd_ea, match_asd_ea_c)
                     )]
scz_asd_ea_gwas_data[, union_match_all := pmax(union_match_scz_ea,
                                               union_match_asd_ea)]
table(scz_asd_ea_gwas_data$union_match_all)
# 1
# 6523119
any(is.na(scz_asd_ea_gwas_data$union_match_all))
# [1] FALSE
# Done no impact from complement matching

# Matches same number as the ASD / SCZ join which is convenient
rm(scz_asd_gwas_data, ea_gwas_data)

# Check the INFO imputation scores for the provided datasets
summary(scz_asd_ea_gwas_data$asd_info)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7000  0.9470  0.9850  0.9564  0.9950  1.0600
summary(scz_asd_ea_gwas_data$scz_info)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.302   0.951   0.986   0.953   0.996   1.050

# Which genomic positions have more than one SNP:
duplicate_hg19_ids <- scz_asd_ea_gwas_data %>%
  .[, n_rs_ids := .N, by = hg19_id] %>%
  .[n_rs_ids > 1, hg19_id] %>%
  unique()
# there are none

# Filter SNPs based on 1000 Genomes MAF -----------------------------------

g1000_eur_fam_path <- "data/g1000_eur/g1000_eur.fam"
g1000_eur_bim_path <- "data/g1000_eur/g1000_eur.bim"
g1000_eur_bed_path <- "data/g1000_eur/g1000_eur.bed"

# Read in the PLINK data using snpStats:
ref_snps_plink <- read.plink(g1000_eur_bed_path,
                             g1000_eur_bim_path,
                             g1000_eur_fam_path)

# Obtain the SnpMatrix object (genotypes) table from ref_snps_plink list
ref_genotypes <- ref_snps_plink$genotypes
print(ref_genotypes)
# A SnpMatrix with  503 rows and  22665064 columns
# Row names:  HG00096 ... NA20832
# Col names:  rs537182016 ... rs781880

#Obtain the SNP information from ref_snps_plink list
ref_genobim <- ref_snps_plink$map
colnames(ref_genobim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
ref_genobim <- as.data.table(ref_genobim)
head(ref_genobim)
#    chr         SNP gen.dist position A1 A2
# 1:   1 rs537182016       NA    10539  A  C
# 2:   1 rs575272151       NA    11008  G  C
# 3:   1 rs544419019       NA    11012  G  C
# 4:   1 rs540538026       NA    13110  A  G
# 5:   1  rs62635286       NA    13116  G  T
# 6:   1 rs200579949       NA    13118  G  A
# Remove raw file to free up memory
rm(ref_snps_plink)

# Create SNP summary statistics (MAF, call rate, etc.)
snp_sum_col <- col.summary(ref_genotypes)

# Set thresholds for MAF:
min_maf <- 0.05
# Filter SNPs on MAF and call rate
use_snps_i <- with(snp_sum_col,
                   (!is.na(MAF) & MAF >= min_maf))
# Remove NA's as well
use_snps_i[is.na(use_snps_i)] <- FALSE
cat(ncol(ref_genotypes) - sum(use_snps_i),
    "SNPs will be removed due to MAF < 0.05.\n")
# 16445832 SNPs will be removed due to MAF < 0.05.
rm(snp_sum_col)

# Subset genotype and SNP summary data for SNPs that pass MAF criteria
ref_genotypes <- ref_genotypes[, use_snps_i]
rm(use_snps_i)

# Now filter ref_genobim
ref_genobim <- ref_genobim[which(ref_genobim$SNP %in% colnames(ref_genotypes)),]
# No longer need the genotype data available:
rm(ref_genotypes)

# Create a column in the reference data information with the unique location
# information for matchign with the GWAS data:
ref_genobim[, hg19_id := paste0("chr", chr, ":", position, "-", position)]
length(unique(ref_genobim$SNP))
# [1] 6219232
length(unique(ref_genobim$hg19_id))
# [1] 6219232 - no duplicates

# Only keep the reference data locations for those that are in the GWAS data:
gwas_ref_genobim <- ref_genobim[hg19_id %in% scz_asd_ea_gwas_data$hg19_id,]
nrow(gwas_ref_genobim)
# [1] 5241136
# Major reduction in the number of SNPs

# Next need to check that these have the same alleles. First rename the allele
# columns in the reference data:
setnames(gwas_ref_genobim, old = "A1", new = "ref_a1")
setnames(gwas_ref_genobim, old = "A2", new = "ref_a2")

# Make the complementary ones:
gwas_ref_genobim[, ':='(ref_ca1 = as.character(allele_complements[ref_a1]),
                        ref_ca2 = as.character(allele_complements[ref_a2]))]
any(is.na(gwas_ref_genobim$ref_ca1))
# [1] FALSE
any(is.na(gwas_ref_genobim$ref_ca2))
# [1] FALSE
# These are all single letter changes

# Only consider the genomic locations in the reference data and check to make
# sure that the alleles match (in either order):
scz_asd_ea_gwas_data <- scz_asd_ea_gwas_data %>%
  .[hg19_id %in% gwas_ref_genobim$hg19_id,] %>%
  merge(gwas_ref_genobim[, .(hg19_id, ref_a1, ref_a2, ref_ca1, ref_ca2)],
        all.x = TRUE, by = "hg19_id") %>%
  # Create new columns denoting if the alleles match (in either direction) for
  # either the observed or complementary:
  .[, ':='(ref_match_o = as.numeric((scz_a1 == ref_a1 & scz_a2 == ref_a2) |
                                      (scz_a1 == ref_a2 & scz_a2 == ref_a1)),
           ref_match_c = as.numeric((scz_a1 == ref_ca1 & scz_a2 == ref_ca2) |
                                      (scz_a1 == ref_ca2 & scz_a2 == ref_ca1)))]
table(scz_asd_ea_gwas_data$ref_match_o)
#       1
# 5241136 - all match amazingly
table(scz_asd_ea_gwas_data$ref_match_c)
# 0       1
# 4449722  791414 - basically a sign of how many are just flipped

# Remove the matching columns since they are no longer necessary:
scz_asd_ea_gwas_data <- scz_asd_ea_gwas_data %>%
  .[, !c("match_scz_asd", "match_scz_asd_c", "union_match_scz_asd",
         "match_scz_ea", "match_scz_ea_c", "match_asd_ea", "match_asd_ea_c",
         "union_match_scz_ea", "union_match_asd_ea", "union_match_all",
         "n_rs_ids", "ref_match_o", "ref_match_c")]

# Save this dataset:
# fwrite(scz_asd_ea_gwas_data,
#        "data/gwas/scz_asd_ea_gwas_results.csv")


# Add hg38 coordinates for joining BrainVar eQTLs -------------------------

# Create a dataset with only the SNPs hg19_id - will then save this
# file as input for liftOver:
hg19_ids <- scz_asd_ea_gwas_data %>%
  dplyr::select(hg19_id) %>%
  distinct()
# Save as liftOver input:
hg19_ids %>%
  write_csv("data/gwas/liftover/input/gwas_snps_hg19_ids.csv",
            col_names = FALSE)

# USE THE LIFTOVER TOOL:

# Next load the errors - to which rows will not have a new ID to assign corresponding to hg19:
snp_errors <- fread("data/gwas/liftover/output/gwas_snps_err.txt",
                    header = FALSE) %>%
  as.data.frame() %>%
  rename(hg19_id = V1) %>%
  filter(str_detect(hg19_id, "chr"))

# Which rows are these:
snp_errors_i <- which(hg19_ids$hg19_id %in% snp_errors$hg19_id)

# Next load the output from liftover:
snps_hg38_ids <- fread("data/gwas/liftover/output/gwas_snps_output.bed",
                       header = FALSE) %>%
  as.data.frame() %>%
  rename(hg38_id = V1)

# Create a new column for hg19_ids hg38_id that is initially all NA:
hg19_ids$hg38_id <- NA
# Now update for the rows that are not errors:
hg19_ids$hg38_id[-snp_errors_i] <- snps_hg38_ids$hg38_id

# Join these data to the original gwas data
scz_asd_ea_gwas_data_merge <- merge(scz_asd_ea_gwas_data, hg19_ids)

# Finally drop the complementary allele columns since they are no longer
# necessary and can just use asd_a1 and asd_a2 since they match across the
# three phenotypes:
scz_asd_ea_gwas_data_merge <- scz_asd_ea_gwas_data_merge %>%
  .[, !c("scz_ca1", "scz_ca2", "asd_ca1", "asd_ca2", "ea_ca1", "ea_ca2",
         "ref_ca1", "ref_ca2", "scz_a1", "scz_a2", "ea_a1", "ea_a2")]

# Save this update:
fwrite(scz_asd_ea_gwas_data_merge,
       "data/gwas/scz_asd_ea_gwas_results.csv")




