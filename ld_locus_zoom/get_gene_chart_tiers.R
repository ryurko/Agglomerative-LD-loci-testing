# PURPOSE: Create helper function for determining gene order based on overlap
#          and the size of the genes with respect to start and end positions

#' Determine the tier for displaying genes
#'
#' @param gene_ids Vector of unique gene identifiers
#' @param start_pos Vector of gene start positions
#' @param end_pos Vector of gene end positions
#' @return Vector of numeric values indicated which tier each gene (with respect)
#' to the input order is displayed on the app visualization.

get_gene_chart_tiers <- function(gene_ids, start_pos, end_pos) {

  # First determine the overlap between genes based on their start / end positions
  gene_table_dt <- data.table::data.table(
    start = start_pos,
    end = end_pos,
    id = gene_ids,
    key = c("start", "end")
  )
  gene_overlap <- data.table::foverlaps(gene_table_dt, gene_table_dt)
  gene_overlap <- dplyr::filter(dplyr::as_tibble(gene_overlap), id != i.id)

  gene_tiers <- rep(0, length(gene_ids))
  names(gene_tiers) <- gene_ids

  for (i in 1:length(gene_ids)) {

    # Get the next gene (these are assumed to be sorted by size)
    gene_i <- gene_ids[i]

    # Does it overlap with any - if not then just assign to 1:
    if (!(gene_i %in% gene_overlap$i.id)) {
      gene_tiers[gene_i] <- 1
    } else {
      # Otherwise need to find the maximum of its overlaps and add 1:
      overlaps <- dplyr::pull(dplyr::filter(gene_overlap, i.id == gene_i),
                              id)
      gene_tiers[gene_i] <- max(gene_tiers[overlaps]) + 1
    }
  }

  return(gene_tiers)
}


