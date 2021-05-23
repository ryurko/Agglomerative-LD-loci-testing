# PURPOSE: Generate the kernel smoothing results for the positional SNPs at
#          both the gene and LD loci level, with null results generated at
#          the LD loci level. Results are generated across all types of
#          SNP-to-gene assignments

# SCZ VERSION

library(tidyverse)

# Source necessary functions
source("code/s7_prep_app_smoothing/kernel_smoothing_fns.R")

# Load the input positional SNP to LD loci and bandwidth datasets --------

# Positional + eSNPS
pos_esnps_loc_snp_ld_loci_results <-
  read_csv("data/kernel_smoothing/scz/input/positional_esnps/positional_snp_ld_loci_data.csv")
pos_esnps_loc_snp_gene_results <-
  read_csv("data/kernel_smoothing/scz/input/positional_esnps/positional_snp_gene_ld_loci_data.csv")
pos_esnps_bw_data <-
  read_csv("data/kernel_smoothing/scz/input/positional_esnps/ld_loci_gcv_bw.csv")

# Positional
pos_loc_snp_ld_loci_results <-
  read_csv("data/kernel_smoothing/scz/input/positional/positional_snp_ld_loci_data.csv")
pos_loc_snp_gene_results <-
  read_csv("data/kernel_smoothing/scz/input/positional/positional_snp_gene_ld_loci_data.csv")
pos_bw_data <-
  read_csv("data/kernel_smoothing/scz/input/positional/ld_loci_gcv_bw.csv")

# Generate locus clusters -------------------------------------------------


# # Size cutoff for Positional + eSNPs LD loci results:
# pos_esnps_loc_snp_ld_loci_results %>%
#   group_by(ld_locus_id) %>%
#   summarize(locus_size = max(bp) - min(bp)) %>%
#   ungroup() %>%
#   pull(locus_size) %>%
#   summary()
# # Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# #    0     10650     70781   4027317    282841 234665644
#
# # Positional results:
# pos_loc_snp_ld_loci_results %>%
#   group_by(ld_locus_id) %>%
#   summarize(locus_size = max(bp) - min(bp)) %>%
#   ungroup() %>%
#   pull(locus_size) %>%
#   summary()
# #    Min.   1st Qu.   Median      Mean    3rd Qu.     Max.
# #        0     12294     72580   3176507    253806 226622871

# Will use the 75th percentile to make the cutoff for clustering
pos_esnps_loc_snp_ld_loci_results <- pos_esnps_loc_snp_ld_loci_results %>%
  group_by(ld_locus_id) %>%
  mutate(intra_locus_cluster = get_intra_loci_groups(bp, 282841)) %>%
  ungroup()
pos_loc_snp_ld_loci_results <- pos_loc_snp_ld_loci_results %>%
  group_by(ld_locus_id) %>%
  mutate(intra_locus_cluster = get_intra_loci_groups(bp, 253806)) %>%
  ungroup()


# Generate the LD loci and gene-level results ----------------------------

# Positional + eSNPs
pos_esnps_ld_loci_level_smoothing <-
  get_locus_cluster_level_smooth_results(pos_esnps_loc_snp_ld_loci_results, pos_esnps_bw_data,
                                   "scz_z_squared")

pos_esnps_gene_level_smoothing <-
  get_gene_level_smooth_results(pos_esnps_loc_snp_gene_results, pos_esnps_bw_data,
                                "scz_z_squared", n_gene_points = 100)
# Join info to gene-level and save:
pos_esnps_gene_level_smoothing <- pos_esnps_gene_level_smoothing %>%
  dplyr::rename(smooth_scz_zsquared = smooth_zsquared) %>%
  dplyr::left_join(pos_esnps_loc_snp_gene_results,
                   by = c("bp", "ld_locus_id", "ensembl_id")) %>%
  mutate(is_fake_bp = as.numeric(is.na(hg38_id))) %>%
  # For ease later on in app, save the snp id with the version's hg coordinates
  dplyr::rename(snp_id = hg38_id)
# And save
write_csv(pos_esnps_gene_level_smoothing,
          "data/kernel_smoothing/scz/output/positional_esnps/positional_gene_level_smoothing.csv")

# Repeat these steps for the remaining assignment types

# Positional
pos_ld_loci_level_smoothing <-
  get_locus_cluster_level_smooth_results(pos_loc_snp_ld_loci_results, pos_bw_data,
                                   "scz_z_squared")

pos_gene_level_smoothing <-
  get_gene_level_smooth_results(pos_loc_snp_gene_results, pos_bw_data,
                                "scz_z_squared", n_gene_points = 100)
pos_gene_level_smoothing <- pos_gene_level_smoothing %>%
  dplyr::rename(smooth_scz_zsquared = smooth_zsquared) %>%
  dplyr::left_join(pos_loc_snp_gene_results,
                   by = c("bp", "ld_locus_id", "ensembl_id")) %>%
  mutate(is_fake_bp = as.numeric(is.na(hg38_id))) %>%
  dplyr::rename(snp_id = hg38_id)
write_csv(pos_gene_level_smoothing,
          "data/kernel_smoothing/scz/output/positional/positional_gene_level_smoothing.csv")


# Generate the null results for each  -------------------------------------

pos_esnps_ld_loci_null_results <-
  get_ld_loci_level_null_results(pos_esnps_loc_snp_ld_loci_results, pos_esnps_bw_data,
                                 n_cores = 4)
pos_ld_loci_null_results <-
  get_ld_loci_level_null_results(pos_loc_snp_ld_loci_results, pos_bw_data,
                                 n_cores = 4)

# Join to the LD loci results and save -----------------------------------

# Join the LD loci level results with the background, null, and LD loci level
# info to then save for use with the LD loci zoom app

# Start with Positional + eSNPs
pos_esnps_ld_loci_level_smoothing <- pos_esnps_ld_loci_level_smoothing %>%
  dplyr::rename(smooth_scz_zsquared = smooth_zsquared) %>%
  # Join the null results
  dplyr::left_join(pos_esnps_ld_loci_null_results,
                   by = c("bp", "ld_locus_id", "intra_locus_cluster")) %>%
  # Join the remaining LD loci info
  dplyr::left_join(pos_esnps_loc_snp_ld_loci_results,
                   by = c("bp", "ld_locus_id", "intra_locus_cluster")) %>%
  dplyr::rename(snp_id = hg38_id) %>%
  mutate(is_fake_bp = as.numeric(is.na(snp_id)))
write_csv(pos_esnps_ld_loci_level_smoothing,
          "data/kernel_smoothing/scz/output/positional_esnps/positional_ld_loci_level_smoothing.csv")

# Repeat for Positional
pos_ld_loci_level_smoothing <- pos_ld_loci_level_smoothing %>%
  dplyr::rename(smooth_scz_zsquared = smooth_zsquared) %>%
  dplyr::left_join(pos_ld_loci_null_results,
                   by = c("bp", "ld_locus_id", "intra_locus_cluster")) %>%
  dplyr::left_join(pos_loc_snp_ld_loci_results,
                   by = c("bp", "ld_locus_id", "intra_locus_cluster")) %>%
  dplyr::rename(snp_id = hg38_id) %>%
  mutate(is_fake_bp = as.numeric(is.na(snp_id)))
write_csv(pos_ld_loci_level_smoothing,
          "data/kernel_smoothing/scz/output/positional/positional_ld_loci_level_smoothing.csv")

