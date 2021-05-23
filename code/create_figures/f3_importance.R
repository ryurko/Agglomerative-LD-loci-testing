# PURPOSE: Create versions of the AdaPT variable importance plots that combines
#          the \pi_1 and \mu results into a single plot, just displaying the
#          top five variables based on both models. Will leave out the AdaPT step
#          explanation in the figures and save that information for the supplement

library(tidyverse)
library(xgboost)
library(cowplot)

# Load AdaPT model data ---------------------------------------------------

pos_esnps_model_data <-
  read_csv("data/adapt_model_data/positional_esnps/rsquared25_model_data.csv")

pos_model_data <-
  read_csv("data/adapt_model_data/positional/rsquared25_model_data.csv")

# Load the AdaPT CV results -----------------------------------------------

# Make a function to do this quicker:
read_xgb_cv_results <- function(assign_type, phenotype, rsquared) {
  readRDS(
    paste0("data/adapt_results/", assign_type, "/",
           phenotype, "/rsquared", rsquared, "_xgb_cv.rds")
  )
}

# Positional + eSNPs

# ASD
pos_esnps_asd_25_xgb_cv <- read_xgb_cv_results("positional_esnps", "asd", 25)

# SCZ
pos_esnps_scz_25_xgb_cv <- read_xgb_cv_results("positional_esnps", "scz", 25)

# EA
pos_esnps_ea_25_xgb_cv <- read_xgb_cv_results("positional_esnps", "ea", 25)

# Positional

# ASD
pos_asd_25_xgb_cv <- read_xgb_cv_results("positional", "asd", 25)

# SCZ
pos_scz_25_xgb_cv <- read_xgb_cv_results("positional", "scz", 25)

# EA
pos_ea_25_xgb_cv <- read_xgb_cv_results("positional", "ea", 25)


# Create datasets summarizing the variable importance  --------------------

# First write a helper function to init the dataset of the variable importance
# across the entire AdaPT search and model types:
create_var_imp_data <- function(adapt_results) {

  map_dfr(1:length(adapt_results$model_fit),
          function(step_i) {

            # Do a try-catch to deal with constants returned by the trees:
            tryCatch({
              xgb.importance(model = adapt_results$model_fit[[step_i]]) %>%
                as_tibble() %>%
                mutate(adapt_step = step_i,
                       model_type = ifelse(step_i %% 2 == 0, "mu", "pi_1"))
            },
            error = function(cond) {
              message(paste0("Did not use any variables for model ", model_type,
                             " at step ", step_i, " (", model_step_i, "-model step)"))
              tibble(Feature = "NONE",
                     adapt_step = step_i,
                     model_type = ifelse(step_i %% 2 == 0, "mu", "pi_1"))
            })
          })
}

# Apply to each set of results, starting with Positional + eSNPs:
pos_esnps_asd_var_imp_data <- create_var_imp_data(pos_esnps_asd_25_xgb_cv)
pos_esnps_scz_var_imp_data <- create_var_imp_data(pos_esnps_scz_25_xgb_cv)
pos_esnps_ea_var_imp_data <- create_var_imp_data(pos_esnps_ea_25_xgb_cv)
# Positional:
pos_asd_var_imp_data <- create_var_imp_data(pos_asd_25_xgb_cv)
pos_scz_var_imp_data <- create_var_imp_data(pos_scz_25_xgb_cv)
pos_ea_var_imp_data <- create_var_imp_data(pos_ea_25_xgb_cv)


# Create summary dataset of all AdaPT steps -------------------------------

# Write a helper pipeline for summarizing the data - assuming phenotype is included:
summarize_adapt_var_imp <- . %>%
  filter(Feature != "none") %>%
  group_by(phenotype) %>%
  mutate(n_model_steps = max(adapt_step) / 2) %>%
  ungroup() %>%
  group_by(phenotype, Feature, model_type) %>%
  summarize(ave_gain = sum(Gain) / first(n_model_steps)) %>%
  ungroup() %>%
  group_by(phenotype, Feature) %>%
  mutate(total_ave_gain = sum(ave_gain)) %>%
  ungroup()

# Generate the Positional + eSNPs summary:
pos_esnps_var_imp_summary <- bind_rows(mutate(pos_esnps_asd_var_imp_data,
                                              phenotype = "asd"),
                                       mutate(pos_esnps_scz_var_imp_data,
                                              phenotype = "scz"),
                                       mutate(pos_esnps_ea_var_imp_data,
                                              phenotype = "ea")) %>%
  summarize_adapt_var_imp()
# Positional:
pos_var_imp_summary <- bind_rows(mutate(pos_asd_var_imp_data,
                                        phenotype = "asd"),
                                 mutate(pos_scz_var_imp_data,
                                        phenotype = "scz"),
                                 mutate(pos_ea_var_imp_data,
                                        phenotype = "ea")) %>%
  summarize_adapt_var_imp()


# Create stacked barcharts for comparison ---------------------------------

# What are the five total variables for each phenotype:
pos_esnps_top_var_imp <- pos_esnps_var_imp_summary %>%
  dplyr::select(phenotype, Feature, total_ave_gain) %>%
  distinct() %>%
  group_by(phenotype) %>%
  top_n(5, wt = total_ave_gain) %>%
  ungroup() %>%
  mutate(is_top = 1)

# Now join to the whole data
pos_esnps_var_imp_plot <- pos_esnps_var_imp_summary %>%
  dplyr::left_join(dplyr::select(pos_esnps_top_var_imp, -total_ave_gain),
                   by = c("phenotype", "Feature")) %>%
  filter(is_top == 1) %>%
  mutate(model_type = ifelse(model_type == "mu", "Non-null effect size",
                             "Non-null probability"),
         phenotype = fct_relevel(toupper(phenotype), "ASD", "SCZ", "EA"),
         Feature = ifelse(Feature == "gtex_any_gene_green",
                          "GTEx WGCNA\ngreen module indicator", Feature),
         Feature = fct_recode(Feature, `EA z statistics` = "cap_ea_quad_z",
                              `SCZ z statistics` = "cap_scz_quad_z",
                              `ASD z statistics` = "cap_asd_quad_z",
                              `Min. LOEUF` = "min_loeuf",
                              `Number of SNPs` = "n_snps",
                              `Ave |pre beta|` = "ave_abs_pre_beta"),
         Feature = tidytext::reorder_within(Feature, total_ave_gain, phenotype)) %>%
  ggplot(aes(x = Feature, y = ave_gain, fill = model_type, color = model_type)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_manual(values = c("darkblue", "darkorange"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = c("darkblue", "darkorange"), guide = FALSE) +
  tidytext::scale_x_reordered() +
  coord_flip() +
  facet_wrap(~ phenotype, ncol = 3, scales = "free") +
  labs(x = "Metadata", y = "Importance",
       fill = "AdaPT model type:") +
  theme_bw() +
  theme(legend.position = c(.15, .2),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14))

save_plot("figures/main/pos_esnps_importance.jpg",
          pos_esnps_var_imp_plot, ncol = 3, nrow = 1)
save_plot("figures/main/pos_esnps_importance.pdf",
          pos_esnps_var_imp_plot, ncol = 3, nrow = 1)

# Repeat for the Positional results
pos_top_var_imp <- pos_var_imp_summary %>%
  dplyr::select(phenotype, Feature, total_ave_gain) %>%
  distinct() %>%
  group_by(phenotype) %>%
  top_n(5, wt = total_ave_gain) %>%
  ungroup() %>%
  mutate(is_top = 1)

# Now join to the whole data
pos_var_imp_plot <- pos_var_imp_summary %>%
  dplyr::left_join(dplyr::select(pos_top_var_imp, -total_ave_gain),
                   by = c("phenotype", "Feature")) %>%
  filter(is_top == 1) %>%
  mutate(model_type = ifelse(model_type == "mu", "Non-null effect size",
                             "Non-null probability"),
         phenotype = fct_relevel(toupper(phenotype), "ASD", "SCZ", "EA"),
         Feature = ifelse(Feature == "gtex_any_gene_green",
                          "GTEx WGCNA\ngreen module indicator", Feature),
         Feature = fct_recode(Feature, `EA z statistics` = "cap_ea_quad_z",
                              `SCZ z statistics` = "cap_scz_quad_z",
                              `ASD z statistics` = "cap_asd_quad_z",
                              `Min. LOEUF` = "min_loeuf",
                              `Number of SNPs` = "n_snps",
                              `Ave |pre beta|` = "ave_abs_pre_beta"),
         Feature = tidytext::reorder_within(Feature, total_ave_gain, phenotype)) %>%
  ggplot(aes(x = Feature, y = ave_gain, fill = model_type, color = model_type)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_manual(values = c("darkblue", "darkorange"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = c("darkblue", "darkorange"), guide = FALSE) +
  tidytext::scale_x_reordered() +
  coord_flip() +
  facet_wrap(~ phenotype, ncol = 3, scales = "free") +
  labs(x = "Metadata", y = "Importance",
       fill = "AdaPT model type:") +
  theme_bw() +
  theme(legend.position = c(.15, .2),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14))

save_plot("figures/suppl/pos_importance.jpg",
          pos_var_imp_plot, ncol = 3, nrow = 1)
save_plot("figures/suppl/pos_importance.pdf",
          pos_var_imp_plot, ncol = 3, nrow = 1)




