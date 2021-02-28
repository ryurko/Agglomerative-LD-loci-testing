# PURPOSE: Examine the top settings using the synthetic SCZ p-values

library(tidyverse)

# Positional + eSNPs ------------------------------------------------------

# r^2 = 0.25
pos_esnp_25_tune_results <-
  read_csv("data/adapt_results/positional_esnps/fake_scz_tune/rsquared25_tune_results.csv")

# What are the best settings here for both levels of depth:
pos_esnp_25_tune_results %>%
  filter(alpha == 0.01) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)

# A tibble: 2 x 10
#     Groups:   max_depth [2]
#     phenotype ld_loci_type                method        variables       alpha n_disc nrounds max_depth   eta gamma
#     <chr>     <chr>                        <chr>            <chr>       <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
#   1 fake_scz  pos_esnps_rsquared25 adapt_xgboost              all        0.01    318     400         1  0.06     0
#   2 fake_scz  pos_esnps_rsquared25 adapt_xgboost              all        0.01    343     200         2  0.05     0

# r^2 = 0.50
pos_esnp_50_tune_results <-
  read_csv("data/adapt_results/positional_esnps/fake_scz_tune/rsquared50_tune_results.csv")

# What are the best settings here for both levels of depth:
pos_esnp_50_tune_results %>%
  filter(alpha == 0.01) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)
# A tibble: 2 x 10
#     Groups:   max_depth [2]
#     phenotype         ld_loci_type        method    variables        alpha n_disc nrounds max_depth   eta gamma
#        <chr>                 <chr>         <chr>         <chr>       <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
#   1 fake_scz  pos_esnps_rsquared50 adapt_xgboost           all        0.01    287     300         1  0.04     0
#   2 fake_scz  pos_esnps_rsquared50 adapt_xgboost           all        0.01    282     150         2  0.06     0



# Positional --------------------------------------------------------------

# r^2 = 0.25
pos_25_tune_results <-
  read_csv("data/adapt_results/positional/fake_scz_tune/rsquared25_tune_results.csv")

# What are the best settings here for both levels of depth:
pos_25_tune_results %>%
  filter(alpha == 0.01) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)

# A tibble: 2 x 10
#     Groups:   max_depth [2]
#     phenotype   ld_loci_type        method     variables      alpha n_disc nrounds max_depth   eta gamma
#         <chr>          <chr>         <chr>         <chr>      <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
#   1 fake_scz  pos_rsquared25 adapt_xgboost          all        0.01    221     100         1  0.02     0
#   2 fake_scz  pos_rsquared25 adapt_xgboost          all        0.01    223     100         2  0.03     0

# r^2 = 0.50
pos_50_tune_results <-
  read_csv("data/adapt_results/positional/fake_scz_tune/rsquared50_tune_results.csv")

# What are the best settings here for both levels of depth:
pos_50_tune_results %>%
  filter(alpha == 0.01) %>%
  group_by(max_depth) %>%
  arrange(desc(n_disc)) %>%
  slice(1)
# A tibble: 2 x 10
#     Groups:   max_depth [2]
#     phenotype   ld_loci_type        method     variables       alpha n_disc nrounds max_depth   eta gamma
#        <chr>           <chr>         <chr>         <chr>       <dbl>  <dbl>   <dbl>     <dbl> <dbl> <dbl>
#   1 fake_scz  pos_rsquared50 adapt_xgboost           all        0.01    221     100         1  0.02     0
#   2 fake_scz  pos_rsquared50 adapt_xgboost           all        0.01    303     100         2  0.02     0

