.blk_apply <- function(iobject, split, v = "cat") {
  o <- match(igraph::vertex.attributes(iobject)$name, split$vertex)
  o_block <- split$block[o]
  blk_return <- igraph::set.vertex.attribute(iobject, v, value = o_block)
  return(blk_return)
}

#' @export
concor_make_igraph <- function(adj_list, nsplit = 1) {
  adj_list <- .concor_validitycheck(adj_list)
  concor_out <- suppressWarnings(concor(adj_list, nsplit))

  all_unweighted <- all(sapply(adj_list, function(x) all(x %in% c(0,1))))
  if (all_unweighted) {
    igraph_list <- lapply(adj_list,
                          function(x) igraph::graph_from_adjacency_matrix(x))
  } else {
    igraph_list <- lapply(adj_list,
                          function(x) igraph::graph_from_adjacency_matrix(x, weighted = "weight"))
  }

  v <- paste("csplit", nsplit, sep = "")
  igraph_out <- lapply(igraph_list, function(x) .blk_apply(x, concor_out, v))

  return(igraph_out)
}

#' @export
concor_igraph_apply <- function(igraph_list, nsplit = 1) {
  any_weighted <- any(sapply(igraph_list, igraph::is_weighted))
  all_weighted <- all(sapply(igraph_list, igraph::is_weighted))

  if (all_weighted) {
    adj_list <- lapply(igraph_list,
                       function(x) igraph::as_adjacency_matrix(x,
                                                               attr = "weight",
                                                               sparse = FALSE))
  } else {
    if (any_weighted) {
      warning("Some but not all input graphs are weighted; ignoring weights.\n")
    }

    adj_list <- lapply(igraph_list,
                       function(x) igraph::as_adjacency_matrix(x,
                                                               sparse = FALSE))
  }

  concor_out <- suppressWarnings(concor(adj_list, nsplit))
  v <- paste("csplit", nsplit, sep = "")
  igraph_out <- lapply(igraph_list, function(x) .blk_apply(x, concor_out, v))

  return(igraph_out)
}

#' @export
plot_socio <- function(iobject, nsplit = NULL) {
  split_name <- paste0("csplit", nsplit)
  vcolors <- igraph::vertex.attributes(iobject)[[split_name]]
  igraph::plot.igraph(iobject, vertex.color = vcolors,
                      vertex.label = NA, vertex.size = 5, edge.arrow.size = .3)
}

