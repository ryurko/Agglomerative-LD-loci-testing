# PURPOSE: Create figure with comparison of AdaPT results to BH and intercept-
#          only results for each of the three phenotypes using the r^2 = 0.75 data

library(tidyverse)
library(cowplot)

# Load AdaPT model data ---------------------------------------------------

pos_esnps_model_data <-
  read_csv("data/adapt_model_data/positional_esnps/rsquared75_model_data.csv")

pos_model_data <-
  read_csv("data/adapt_model_data/positional/rsquared75_model_data.csv")

# Load the intercept-only AdaPT results -----------------------------------

# Make a function to do this quicker:
read_int_only_results <- function(assign_type, phenotype, rsquared) {
  readRDS(
    paste0("data/adapt_results/", assign_type, "/",
           phenotype, "_rsquared", rsquared, "_int_only.rds")
  )
}

# Positional + eSNPs

# ASD
pos_esnps_asd_75_int_only <- read_int_only_results("positional_esnps", "asd", 75)

# SCZ
pos_esnps_scz_75_int_only <- read_int_only_results("positional_esnps", "scz", 75)

# EA
pos_esnps_ea_75_int_only <- read_int_only_results("positional_esnps", "ea", 75)

# Positional

# ASD
pos_asd_75_int_only <- read_int_only_results("positional", "asd", 75)

# SCZ
pos_scz_75_int_only <- read_int_only_results("positional", "scz", 75)

# EA
pos_ea_75_int_only <- read_int_only_results("positional", "ea", 75)


# Load the AdaPT CV results -----------------------------------------------

# Make a function to do this quicker:
read_xgb_cv_results <- function(assign_type, phenotype, rsquared) {
  readRDS(
    paste0("data/adapt_results/", assign_type, "/",
           phenotype, "_rsquared", rsquared, "_xgb_cv.rds")
  )
}

# Positional + eSNPs

# ASD
pos_esnps_asd_75_xgb_cv <- read_xgb_cv_results("positional_esnps", "asd", 75)

# SCZ
pos_esnps_scz_75_xgb_cv <- read_xgb_cv_results("positional_esnps", "scz", 75)

# EA
pos_esnps_ea_75_xgb_cv <- read_xgb_cv_results("positional_esnps", "ea", 75)

# Positional

# ASD
pos_asd_75_xgb_cv <- read_xgb_cv_results("positional", "asd", 75)

# SCZ
pos_scz_75_xgb_cv <- read_xgb_cv_results("positional", "scz", 75)

# EA
pos_ea_75_xgb_cv <- read_xgb_cv_results("positional", "ea", 75)


# Create dot-plots comparing number of discoveries by each ----------------

# Make a function that returns the number of BH discoveries for each vector:
return_bh_n_disc <- function(model_data_list, pval_type, target_alpha = 0.05) {
  map_dbl(model_data_list,
          function(model_data) {
            sum(p.adjust(model_data[[pval_type]],
                         method = "BH") <= target_alpha)
          })
}

return_adapt_n_disc <- function(adapt_model_list, target_alpha = 0.05) {
  map_dbl(adapt_model_list,
          function(adapt_model) {
            sum(adapt_model$qvals <= target_alpha)
          })
}

# Create the dotplot comparing the number of discoveries

adapt_dotplot <- bind_rows({
  tibble(type = c("Positional + eSNPs", "Positional"),
         BH = return_bh_n_disc(list(pos_esnps_model_data,
                                    pos_model_data),
                               pval_type = "asd_vegas_pval"),
         `AdaPT: intercept-only` =
           return_adapt_n_disc(
             list(pos_esnps_asd_75_int_only,
                  pos_asd_75_int_only)),
         `AdaPT: XGBoost` =
           return_adapt_n_disc(
             list(pos_esnps_asd_75_xgb_cv,
                  pos_asd_75_xgb_cv))) %>%
    pivot_longer(BH:`AdaPT: XGBoost`,
                 names_to = "method",
                 values_to = "n_disc") %>%
    mutate(method = fct_relevel(method, "BH"),
           type = fct_relevel(type, "Positional + eSNPs", "Positional"),
           phenotype = "ASD")
}, {
  tibble(type = c("Positional + eSNPs", "Positional"),
         BH = return_bh_n_disc(list(pos_esnps_model_data,
                                    pos_model_data),
                               pval_type = "scz_vegas_pval"),
         `AdaPT: intercept-only` =
           return_adapt_n_disc(
             list(pos_esnps_scz_75_int_only,
                  pos_scz_75_int_only)),
         `AdaPT: XGBoost` =
           return_adapt_n_disc(
             list(pos_esnps_scz_75_xgb_cv,
                  pos_scz_75_xgb_cv))) %>%
    pivot_longer(BH:`AdaPT: XGBoost`,
                 names_to = "method",
                 values_to = "n_disc") %>%
    mutate(method = fct_relevel(method, "BH"),
           type = fct_relevel(type, "Positional + eSNPs", "Positional"),
           phenotype = "SCZ")
},
{
  tibble(type = c("Positional + eSNPs", "Positional"),
         BH = return_bh_n_disc(list(pos_esnps_model_data,
                                    pos_model_data),
                               pval_type = "ea_vegas_pval"),
         `AdaPT: intercept-only` =
           return_adapt_n_disc(
             list(pos_esnps_ea_75_int_only,
                  pos_ea_75_int_only)),
         `AdaPT: XGBoost` =
           return_adapt_n_disc(
             list(pos_esnps_ea_75_xgb_cv,
                  pos_ea_75_xgb_cv))) %>%
    pivot_longer(BH:`AdaPT: XGBoost`,
                 names_to = "method",
                 values_to = "n_disc") %>%
    mutate(method = fct_relevel(method, "BH"),
           type = fct_relevel(type, "Positional + eSNPs", "Positional"),
           phenotype = "EA")
}) %>%
  mutate(phenotype = fct_relevel(phenotype, "ASD", "SCZ", "EA"),
         type = fct_rev(type)) %>%
  ggplot(aes(x = type, y = n_disc, fill = method)) +
  geom_bar(width = .2, stat = "identity",
           position = position_dodge(width = .75)) +
  geom_point(position = position_dodge(width = .75),
             size = 10, shape = 21, color = "white") +
  theme_bw() +
  ggthemes::scale_fill_colorblind() +
  ggthemes::scale_color_colorblind(guide = FALSE) +
  labs(x = "SNP-gene assignment", y = "Number of detected genes/loci",
       fill = "Method") +
  coord_flip() +
  facet_wrap(~ phenotype, ncol = 3, scales = "free_x") +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 26),
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        plot.margin = margin(2, 2, 2, 2, "cm"))


# Save --------------------------------------------------------------------

save_plot("figures/suppl/f_rsquared75_adapt_disc_plot.jpg",
          adapt_dotplot, ncol = 3, nrow = 2)
save_plot("figures/suppl/f_rsquared75_adapt_disc_plot.pdf",
          adapt_dotplot, ncol = 3, nrow = 2)

