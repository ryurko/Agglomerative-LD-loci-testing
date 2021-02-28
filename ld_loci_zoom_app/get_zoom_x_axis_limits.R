# PURPOSE: Create helper function for returning the plot axes based on the
#          loci smoothing and gene positions

#' Return axes for zoom display
#'
#' @param loci_smoothing_bp Vector of BP positions for the loci-level
#' smoothing data
#' @param gene_starts Vector of gene start positions
#' @param gene_ends Vector of gene end positions
#' @return List with chart_x_min and chart_x_max positions to use for both the
#' smoothing and gene display

get_zoom_x_axis_limits <- function(loci_smoothing_bp, gene_starts, gene_ends) {
  bp_min <- min(loci_smoothing_bp)
  bp_max <- max(loci_smoothing_bp)
  gene_min <- min(gene_starts)
  gene_max <- max(gene_ends)
  chart_x_min <- pmin(bp_min, gene_min)
  chart_x_max <- pmax(bp_max, gene_max)
  list("chart_x_min" = chart_x_min, "chart_x_max" = chart_x_max)
}
