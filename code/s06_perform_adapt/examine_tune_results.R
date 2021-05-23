# PURPOSE: Examine the tuning results

library(tidyverse)

# Helper function to paste for results:
get_tune_file_path <- function(assign_type, ld_threshold) {
  paste0("data/adapt_results/", assign_type,
         "/fake_scz_tune/rsquared", ld_threshold * 100,
         "_tune_results.csv")
}


# Positional + eSNPs ------------------------------------------------------


# LD threshold = 0.25 -----------------------------------------------------

read_csv(get_tune_file_path("positional_esnps", 0.25)) %>%
  filter(alpha == 0.05) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)

# # A tibble: 2 x 10
# # Groups:   max_depth [2]
#     phenotype ld_loci_type         method        variables alpha n_disc nrounds max_depth   eta gamma
#     <chr>     <chr>                <chr>         <chr>     <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
# 1   fake_scz  pos_esnps_rsquared25 adapt_xgboost all        0.05    729     400         1  0.06     0
# 2   fake_scz  pos_esnps_rsquared25 adapt_xgboost all        0.05    711     250         2  0.05     0

# LD threshold = 0.50 -----------------------------------------------------

read_csv(get_tune_file_path("positional_esnps", 0.5)) %>%
  filter(alpha == 0.05) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)

# # A tibble: 2 x 10
# # Groups:   max_depth [2]
#     phenotype ld_loci_type         method        variables alpha n_disc nrounds max_depth   eta gamma
#     <chr>     <chr>                <chr>         <chr>     <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
#   1 fake_scz  pos_esnps_rsquared50 adapt_xgboost all        0.05    763     450         1  0.06     0
#   2 fake_scz  pos_esnps_rsquared50 adapt_xgboost all        0.05    784     300         2  0.05     0

# LD threshold = 0.75 -----------------------------------------------------


read_csv(get_tune_file_path("positional_esnps", 0.75)) %>%
  filter(alpha == 0.05) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)
# # A tibble: 2 x 10
# # Groups:   max_depth [2]
#     phenotype ld_loci_type         method        variables alpha n_disc nrounds max_depth   eta gamma
#     <chr>     <chr>                <chr>         <chr>     <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
#   1 fake_scz  pos_esnps_rsquared75 adapt_xgboost all        0.05    840     450         1  0.05     0
#   2 fake_scz  pos_esnps_rsquared75 adapt_xgboost all        0.05    815     250         2  0.04     0

# Positional --------------------------------------------------------------


# LD threshold = 0.25 -----------------------------------------------------

read_csv(get_tune_file_path("positional", 0.25)) %>%
  filter(alpha == 0.05) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)
# A tibble: 2 x 10
# Groups:   max_depth [2]
#     phenotype ld_loci_type   method        variables alpha n_disc nrounds max_depth   eta gamma
#     <chr>     <chr>          <chr>         <chr>     <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
#   1 fake_scz  pos_rsquared25 adapt_xgboost all        0.05    704     450         1  0.06     0
#   2 fake_scz  pos_rsquared25 adapt_xgboost all        0.05    721     100         2  0.04     0

# LD threshold = 0.50 -----------------------------------------------------


read_csv(get_tune_file_path("positional", 0.5)) %>%
  filter(alpha == 0.05) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)

# A tibble: 2 x 10
# Groups:   max_depth [2]
#     phenotype ld_loci_type   method        variables alpha n_disc nrounds max_depth   eta gamma
#     <chr>     <chr>          <chr>         <chr>     <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
#   1 fake_scz  pos_rsquared50 adapt_xgboost all        0.05    662     150         1  0.06     0
#   2 fake_scz  pos_rsquared50 adapt_xgboost all        0.05    753     250         2  0.06     0

# LD threshold = 0.75 -----------------------------------------------------


read_csv(get_tune_file_path("positional", 0.75)) %>%
  filter(alpha == 0.05) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)

# # A tibble: 2 x 10
# # Groups:   max_depth [2]
#     phenotype ld_loci_type   method        variables alpha n_disc nrounds max_depth   eta gamma
#     <chr>     <chr>          <chr>         <chr>     <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
#   1 fake_scz  pos_rsquared75 adapt_xgboost all        0.05    813     450         1  0.05     0
#   2 fake_scz  pos_rsquared75 adapt_xgboost all        0.05    812     200         2  0.05     0

