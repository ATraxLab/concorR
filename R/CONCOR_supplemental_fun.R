.blk_apply <- function(iobject, split, v = "cat") {
  o <- match(igraph::vertex.attributes(iobject)$name, split$vertex)
  o_block <- split$block[o]
  blk_return <- igraph::set.vertex.attribute(iobject, v, value = o_block)
  return(blk_return)
}

#' @export
concor_make_igraph <- function(adj_list, nsplit = 1) {
  concor_out <- suppressWarnings(concor(adj_list, nsplit))

  igraph_list <- lapply(adj_list,
                        function(x) igraph::graph_from_adjacency_matrix(x))
  v <- paste("csplit", nsplit, sep = "")
  igraph_out <- lapply(igraph_list, function(x) .blk_apply(x, concor_out, v))

  return(igraph_out)
}

#' @export
concor_igraph_apply <- function(igraph_list, nsplit = 1) {
  adj_list <- lapply(igraph_list,
                     function(x) igraph::get.adjacency(x, sparse = FALSE))

  concor_out <- suppressWarnings(concor(adj_list, nsplit))
  v <- paste("csplit", nsplit, sep = "")
  igraph_out <- lapply(igraph_list, function(x) .blk_apply(x, concor_out, v))

  return(igraph_out)
}

#' @export
plot_socio <- function(iobject, nsplit = NULL) {
  split_name <- paste0("csplit", nsplit)
  if(is.null(nsplit)){
    vcolors <- igraph::vertex.attributes(iobject)[[split_name]]
  }else{
    colpal <- viridis::viridis(max(igraph::vertex.attributes(iobject)[[split_name]]))
    vcolors <- colpal[igraph::vertex.attributes(iobject)[[split_name]]]
  }
  igraph::plot.igraph(iobject, vertex.color = vcolors,
                      vertex.label = NA, vertex.size = 5, edge.arrow.size = .3)
}

# library(concorR)
# library(igraph)
# a <- matrix(c(0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0), ncol = 5)
# rownames(a) <- letters[1:5]
# colnames(a) <- letters[1:5]
# g = concor_make_igraph(list(a), nsplit=1)
# gsimp = set_vertex_attr(g[[1]], "csplit1", value = V(g[[1]])$csplit1)
# plot_socio(gsimp,nsplit=1)
