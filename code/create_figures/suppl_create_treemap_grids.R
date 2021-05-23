# PURPOSE: Arrange the SCZ and EA treemap figures into a grid

library(tidyverse)
library(cowplot)

# Get the figure file paths -----------------------------------------------

# SCZ:
scz_bp_match_path <- "figures/suppl/f_scz_bp_revigo_treemap.pdf"
scz_cc_match_path <- "figures/suppl/f_scz_cc_revigo_treemap.pdf"
scz_mf_match_path <- "figures/suppl/f_scz_mf_revigo_treemap.pdf"

# EA:
ea_bp_match_path <- "figures/suppl/f_ea_bp_revigo_treemap.pdf"
ea_cc_match_path <- "figures/suppl/f_ea_cc_revigo_treemap.pdf"
ea_mf_match_path <- "figures/suppl/f_ea_mf_revigo_treemap.pdf"


# Make the arrangements ---------------------------------------------------

# SCZ -----

scz_grid <-
  plot_grid(ggdraw() + draw_image(scz_bp_match_path, scale = 0.95),
            ggdraw() + draw_image(scz_cc_match_path, scale = 0.95),
            ggdraw() + draw_image(scz_mf_match_path, scale = 0.95),
            labels = c("A", "B", "C"), label_fontface = "plain",
            hjust = -0.75, vjust = 1.25, label_size = 24,
            nrow = 3)
save_plot("figures/suppl/f_scz_treemap_grid.jpg",
          scz_grid, ncol = 1, nrow = 3)
save_plot("figures/suppl/f_scz_treemap_grid.pdf",
          scz_grid, ncol = 1, nrow = 3)

# EA -----

ea_grid <-
  plot_grid(ggdraw() + draw_image(ea_bp_match_path, scale = 0.95),
            ggdraw() + draw_image(ea_cc_match_path, scale = 0.95),
            ggdraw() + draw_image(ea_mf_match_path, scale = 0.95),
            labels = c("A", "B", "C"), label_fontface = "plain",
            hjust = -0.75, vjust = 1.25, label_size = 24,
            nrow = 3)
save_plot("figures/suppl/f_ea_treemap_grid.jpg",
          ea_grid, ncol = 1, nrow = 3)
save_plot("figures/suppl/f_ea_treemap_grid.pdf",
          ea_grid, ncol = 1, nrow = 3)


