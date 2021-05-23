# PURPOSE: Combine the various gene-set zoom plots into single figure displays

library(tidyverse)
library(cowplot)

# Make the highlight example ----------------------------------------------

# Get each of their file paths separately
chr17_inv_a_path <- "figures/suppl/f_chr17_inv_reg_all.png"
chr17_inv_b_path <- "figures/suppl/f_chr17_inv_reg_sub1.png"
chr17_inv_c_path <- "figures/suppl/f_chr17_inv_reg_sub2.png"

# Make the arrangement:
chr17_inv_grid <-
  plot_grid(ggdraw() + draw_image(chr17_inv_a_path, scale = 0.95),
            ggdraw() + draw_image(chr17_inv_b_path, scale = 0.95),
            ggdraw() + draw_image(chr17_inv_c_path, scale = 0.95),
            labels = c("A", "B", "C"), label_fontface = "plain",
            hjust = -0.5, vjust = 1.25, label_size = 24,
            nrow = 3)
save_plot("figures/suppl/f_chr17_inv_zoom.jpg",
          chr17_inv_grid, ncol = 1, nrow = 3, base_asp = 3.5)
save_plot("figures/suppl/f_chr17_inv_zoom.pdf",
          chr17_inv_grid, ncol = 1, nrow = 3, base_asp = 3.5)


