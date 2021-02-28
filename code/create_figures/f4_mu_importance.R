# PURPOSE: Create variable importance and partial dependence plots

library(tidyverse)
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
           phenotype, "_rsquared", rsquared, "_xgb_cv.rds")
  )
}

# Starting with r^2 = 0.25

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

# Next make variable importance plots -------------------------------------

library(xgboost)
library(pdp)
library(cowplot)

# Write functions for creating the variable importance and pdp data:
create_var_imp_data <- function(adapt_results, model_type = 1) {
  model_i <- seq(model_type, length(adapt_results$model_fit), by = 2)
  
  map_dfr(1:length(model_i),
          function(step_i) {
            model_step_i <- model_i[step_i]
            
            # Do a try-catch to deal with constants returned by the trees:
            tryCatch({
              xgb.importance(model = adapt_results$model_fit[[model_step_i]]) %>%
                as_tibble() %>%
                mutate(adapt_step = step_i)
            },
            error = function(cond) {
              message(paste0("Did not use any variables for model ", model_type,
                             " at step ", step_i, " (", model_step_i, "-model step)"))
              tibble(Feature = "NONE",
                     adapt_step = step_i)
            })
          })
}


# Define a function to create the variable importance plot for top n variables
# averaged over the steps in the search (just assumed to have all of the correct
# columns):
create_var_imp_plot <- function(var_imp_data, top_n_vars = 5) {
  var_imp_data %>%
    filter(Feature != "none") %>%
    group_by(Feature) %>%
    mutate(ave_gain = sum(Gain, na.rm = FALSE) / max(var_imp_data$adapt_step)) %>%
    ungroup() %>%
    mutate(ave_step_gain = Gain / max(var_imp_data$adapt_step)) %>%
    filter(Feature %in%
             {
               var_imp_data %>%
                 filter(Feature != "none") %>%
                 group_by(Feature) %>%
                 summarize(ave_gain = sum(Gain, na.rm = FALSE) / max(var_imp_data$adapt_step)) %>%
                 dplyr::arrange(desc(ave_gain)) %>%
                 dplyr::slice(1:top_n_vars) %>%
                 pull(Feature)
             }) %>%
    mutate(Feature = fct_reorder(Feature, ave_gain)) %>%
    ggplot(aes(x = Feature, y = ave_step_gain, fill = adapt_step)) +
    geom_bar(stat = "identity", color = "white", size = 0.1,
             width = .5) +
    scale_fill_gradient(low = "darkblue", high = "darkorange",
                        breaks = c(1, max(var_imp_data$adapt_step)),
                        labels = c("Start", "End")) +
    theme_bw() +
    coord_flip() +
    labs(x = "Metadata", y = "Average importance",
         fill = "AdaPT model fitting iteration") +
    theme(legend.position = "bottom")
}

# Function to return list of the plots:

return_adapt_var_imp_plots <- function(adapt_model_list, model_type, top_n_vars = 5) {
  map(adapt_model_list,
      function(adapt_model) {
        model_data <- create_var_imp_data(adapt_model, model_type)
        create_var_imp_plot(model_data, top_n_vars)
      })
}

# Plot variable importance function:
plot_var_imp_grid <- function(var_imp_list) {
  plot_grid(
    plot_grid(plotlist = lapply(1:length(var_imp_list),
                                function(plot_i) {
                                  var_imp_list[[plot_i]] +
                                    theme(legend.position = "none")
                                }),
              ncol = 2, align = "hv", hjust = 0, vjust = 1.25,
              labels = c("Positional", "Positional + eSNPs")),
    get_legend(var_imp_list[[1]]), ncol = 1,
    rel_heights = c(1, .1))
}

# ASD

asd_pi_imp_plots <- return_adapt_var_imp_plots(
  list(pos_asd_25_xgb_cv, pos_esnps_asd_25_xgb_cv), model_type = 1,
  top_n_vars = 5)
plot_var_imp_grid(asd_pi_imp_plots)

asd_mu_imp_plots <- return_adapt_var_imp_plots(
  list(pos_asd_25_xgb_cv, pos_esnps_asd_25_xgb_cv), model_type = 2,
  top_n_vars = 5)
plot_var_imp_grid(asd_mu_imp_plots)

# SCZ

scz_pi_imp_plots <- return_adapt_var_imp_plots(
  list(pos_scz_25_xgb_cv, pos_esnps_scz_25_xgb_cv), model_type = 1,
  top_n_vars = 5)
plot_var_imp_grid(scz_pi_imp_plots)

scz_mu_imp_plots <- return_adapt_var_imp_plots(
  list(pos_scz_25_xgb_cv, pos_esnps_scz_25_xgb_cv), model_type = 2,
  top_n_vars = 5)
plot_var_imp_grid(scz_mu_imp_plots)

# EA

ea_pi_imp_plots <- return_adapt_var_imp_plots(
  list(pos_ea_25_xgb_cv, pos_esnps_ea_25_xgb_cv), model_type = 1,
  top_n_vars = 5)
plot_var_imp_grid(ea_pi_imp_plots)

ea_mu_imp_plots <- return_adapt_var_imp_plots(
  list(pos_ea_25_xgb_cv, pos_esnps_ea_25_xgb_cv), model_type = 2,
  top_n_vars = 5)
plot_var_imp_grid(ea_mu_imp_plots)


# Create a display of the top five variables for each phenotype's \mu models

mu_imp_plot_grid <-
  #plot_grid(
  plot_grid(asd_mu_imp_plots[[2]] +
              labs(x = "Metadata") +
              scale_x_discrete(labels = rev(c("EA z statistics", "SCZ z statistics",
                                              "Number of SNPs", "Min. LOEUF",
                                              "Ave |pre beta|"))) +
              theme(legend.position = "none"),
            scz_mu_imp_plots[[2]] +
              labs(x = "Metadata", fill = "AdaPT search") +
              guides(fill = guide_colourbar(title.vjust = 1)) +
              scale_x_discrete(labels = rev(c("EA z statistics", "Min. LOEUF",
                                              "Number of SNPs", "ASD z statistics",
                                              "GTEx WGCNA\ngreen module indicator"))) +
              theme(legend.position = c(.6, .2),
                    legend.direction = "horizontal",
                    axis.title.y = element_blank()),
            ea_mu_imp_plots[[2]] +
              labs(x = "Metadata") +
              scale_x_discrete(labels = rev(c("Number of SNPs", "Min. LOEUF",
                                              "SCZ z statistics",
                                              "ASD z statistics",
                                              "GTEx WGCNA\ngreen module indicator"))) +
              theme(legend.position = "none", axis.title.y = element_blank()),
            ncol = 3, labels = c("ASD", "SCZ", "EA"),  hjust = 0, vjust = 1.25
  )


save_plot("figures/main/f4_mu_importance.jpg",
          mu_imp_plot_grid, ncol = 3, nrow = 1)
save_plot("figures/main/f4_mu_importance.pdf",
          mu_imp_plot_grid, ncol = 3, nrow = 1)



# Display partial dependence plot examples --------------------------------
library(latex2exp)
create_var_pdp_data <- function(adapt_results, model_type = 1,
                                model_data, model_var) {
  model_i <- seq(model_type, length(adapt_results$model_fit), by = 2)
  
  map_dfr(1:length(model_i),
          function(step_i) {
            model_step_i <- model_i[step_i]
            suppressMessages(pdp::partial(adapt_results$model_fit[[model_step_i]],
                                          pred.var = model_var,
                                          ice = FALSE, prob = as.numeric(model_type == 1),
                                          center = FALSE, plot = FALSE,
                                          quantiles = TRUE, probs = seq(0, 1, by = .025),
                                          train = data.matrix(
                                            model_data[,adapt_results$model_fit[[1]]$feature_names]))) %>%
              as_tibble() %>%
              mutate(adapt_step = step_i)
          })
}


create_var_pdp_plot <- function(var_pdp_data, model_data, pdp_var, type = 1) {
  var_pdp_plot <- var_pdp_data %>%
    ggplot(aes_string(pdp_var)) +
    geom_line(aes(y = yhat,
                  color = adapt_step,
                  group = as.factor(adapt_step))) +
    cowplot::theme_cowplot() +
    scale_color_gradient(low = "darkblue", high = "darkorange",
                         breaks = c(1, max(var_pdp_data$adapt_step)),
                         labels = c("Start", "End")) +
    geom_rug(data = model_data,
             y = rep(1, nrow(model_data)),
             sides = "b",
             color = "black", alpha = 0.15) +
    guides(color = guide_colourbar(barwidth = 10, barheight = .5)) +
    labs(x = pdp_var,
         y = TeX('$\\mu$'),
         color = "AdaPT model fitting iteration") +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8))
  
  if (type == 1) {
    var_pdp_plot <- var_pdp_plot +
      scale_y_continuous(limits = c(0, 1)) +
      labs(y = TeX('$\\pi_1$'))
  }
  var_pdp_plot
}

return_adapt_var_pdp_plots <- function(adapt_model_list,
                                       adapt_data_list,
                                       model_type,
                                       model_var) {
  map(1:length(adapt_model_list),
      function(adapt_i) {
        pdp_data <- create_var_pdp_data(adapt_model_list[[adapt_i]],
                                        model_type,
                                        adapt_data_list[[adapt_i]],
                                        model_var)
        create_var_pdp_plot(pdp_data, adapt_data_list[[adapt_i]],
                            model_var, model_type)
      })
}


# ASD

# EA z-stats:

asd_scz_z_stat_pi_plots <-
  return_adapt_var_pdp_plots(
    list(pos_asd_25_xgb_cv, pos_esnps_asd_25_xgb_cv),
    list(pos_model_data, pos_esnps_model_data),
    1,
    "cap_scz_quad_z")
asd_scz_z_stat_pi_plots[[1]]
asd_scz_z_stat_pi_plots[[2]]

asd_ea_z_stat_mu_plots <-
  return_adapt_var_pdp_plots(
    list(pos_asd_25_xgb_cv, pos_esnps_asd_25_xgb_cv),
    list(pos_model_data, pos_esnps_model_data),
    2,
    "cap_ea_quad_z")
asd_ea_z_stat_pi_plots[[1]]
asd_ea_z_stat_pi_plots[[2]]

asd_ea_z_stat_mu_plots <-
  return_adapt_var_pdp_plots(
    list(pos_asd_25_xgb_cv, pos_esnps_asd_25_xgb_cv),
    list(pos_model_data, pos_esnps_model_data),
    2,
    "cap_ea_quad_z")
asd_ea_z_stat_mu_plots[[1]]
asd_ea_z_stat_mu_plots[[2]]



# Number of SNPs
asd_n_snps_stat_pi_plots <-
  return_adapt_var_pdp_plots(
    list(pos_asd_25_xgb_cv, pos_esnps_asd_25_xgb_cv),
    list(pos_model_data, pos_esnps_model_data),
    1,
    "n_snps")
asd_n_snps_stat_pi_plots[[1]]
asd_n_snps_stat_pi_plots[[2]]

asd_n_snps_stat_mu_plots <-
  return_adapt_var_pdp_plots(
    list(pos_asd_25_xgb_cv, pos_esnps_asd_25_xgb_cv),
    list(pos_model_data, pos_esnps_model_data),
    2,
    "n_snps")
asd_n_snps_stat_mu_plots[[1]]
asd_n_snps_stat_mu_plots[[2]]




# LOEUF
asd_loeuf_stat_pi_plots <-
  return_adapt_var_pdp_plots(
    list(pos_asd_25_xgb_cv, pos_esnps_asd_25_xgb_cv),
    list(pos_model_data, pos_esnps_model_data),
    1,
    "min_loeuf")
asd_loeuf_stat_pi_plots[[1]]
asd_loeuf_stat_pi_plots[[2]]

asd_loeuf_stat_mu_plots <-
  return_adapt_var_pdp_plots(
    list(pos_asd_25_xgb_cv, pos_esnps_asd_25_xgb_cv),
    list(pos_model_data, pos_esnps_model_data),
    2,
    "min_loeuf")
asd_loeuf_stat_mu_plots[[1]]
asd_loeuf_stat_mu_plots[[2]]


# Make a display with the partial dependence plots for SCZ z-stats and number of SNPs
var_pdp_plot_grid <- plot_grid(asd_scz_z_stat_pi_plots[[2]] +
                                 labs(x = "SCZ z statistics") +
                                 theme(legend.position = "none"),
                               asd_n_snps_stat_mu_plots[[2]] +
                                 labs(x = "Number of SNPs",
                                      color = "AdaPT search") +
                                 guides(colour = guide_colourbar(title.vjust = 1, barwidth = 10)) +
                                 theme(legend.position = c(.3, .2),
                                       legend.direction = "horizontal",
                                       legend.text = element_text(size = 8)),
                               ncol = 2, labels = c("A", "B"),  hjust = -0.5, vjust = 1.25
)

save_plot("figures/suppl/f_example_pdp.jpg",
          var_pdp_plot_grid, ncol = 2, nrow = 1)
save_plot("figures/suppl/f_example_pdp.pdf",
          var_pdp_plot_grid, ncol = 2, nrow = 1)


