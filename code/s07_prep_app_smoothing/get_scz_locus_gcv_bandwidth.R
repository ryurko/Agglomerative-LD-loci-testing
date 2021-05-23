# PURPOSE: Determine the GCV bandwidths for the LD loci for each set of results

# SCZ VERSION

library(tidyverse)
library(np)

# Load the positional SNP to LD loci tables ------------------------------

# Positional + eSNPS
pos_esnps_loc_snp_ld_loci_results <-
  read_csv("data/kernel_smoothing/scz/input/positional_esnps/positional_snp_ld_loci_data.csv")
# Positional
pos_loc_snp_ld_loci_results <-
  read_csv("data/kernel_smoothing/scz/input/positional/positional_snp_ld_loci_data.csv")

# Find GCV bandwidths for each LD loci -----------------------------------

# Write a helper function to apply to each dataset for returning a table of
# GCV bandwidths to use for each LD loci:
get_gcv_bw <- function(snp_ld_loci_results) {

  # Assume input dataset has the ld_locus_id, scz_z_squared, and bp columns:
  map_dfr(unique(snp_ld_loci_results$ld_locus_id),
          function(ld_loci_i) {

            ld_loci_data <- snp_ld_loci_results %>%
              filter(ld_locus_id == ld_loci_i)

            # If the LD loci has more than 2 SNPs proceed with GCV to
            # pick bandwidth using all SNPs:
            if (nrow(ld_loci_data) > 2) {
              ld_loci_bw0 <- npregbw(scz_z_squared ~ bp, data = ld_loci_data)$bw

            } else {
              # Otherwise its NA:
              ld_loci_bw0 <- NA
            }

            # Return a tibble with these bandwidths:
            tibble(ld_locus_id = ld_loci_i,
                   gcv_bw = ld_loci_bw0)

          })

}


# Now apply to each dataset:

# Positional + eSNPS
pos_esnps_loc_gcv <- get_gcv_bw(pos_esnps_loc_snp_ld_loci_results)

# Positional
pos_loc_gcv <- get_gcv_bw(pos_loc_snp_ld_loci_results)


# And save each:
write_csv(pos_esnps_loc_gcv,
          "data/kernel_smoothing/scz/input/positional_esnps/ld_loci_gcv_bw.csv")
write_csv(pos_loc_gcv,
          "data/kernel_smoothing/scz/input/positional/ld_loci_gcv_bw.csv")

