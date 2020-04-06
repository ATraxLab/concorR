#CONCOR supplementary functions
#Tyme Suda


.name <- function(mat) {
  a <- 1:nrow(mat)
  vnames <- sprintf("v%03d", a)
  colnames(mat) <- vnames
  rownames(mat) <- vnames
  return(mat)
}

concor_validitycheck <- function(m_list) {
  a <- m_list[[1]]
  for (i in 1:length(m_list)) {
    if (length(a) != length(m_list[[i]])) {
      stop('Adjacency matrixes of mismatched sizes')
    }
  }

  b <- sapply(m_list, function(x) is.null(colnames(x)))
  if (all(b)) {
    warning("node names don't exist\nAdding default node names\n")
    m_list <- lapply(m_list, function(x) .name(x))
    b <- sapply(m_list, function(x) is.null(colnames(x)))
  }
  if (any(b)) {
    stop("Node name mismatch")
  }

  a <- m_list[[1]]
  for (i in 1:length(m_list)) {
    if (all(colnames(a) != colnames(m_list[[i]]))) {
      stop("Node name mismatch")
    }
  }
  return(m_list)
}

blk_apply <- function(iobject, split, v = "cat") {
  o <- match(igraph::vertex.attributes(iobject)$name, split$vertex)
  o_block <- split$block[o]
  blk_return <- igraph::set.vertex.attribute(iobject, v, value = o_block)
  return(blk_return)
}

make_igraph <- function(adj_list, nsplit = 1) {
  concor_out <- suppressWarnings(concor(adj_list, p = nsplit))

  igraph_list <- lapply(adj_list, function(x) igraph::graph_from_adjacency_matrix(x))
  v <- paste("csplit", nsplit, sep = "")
  igraph_out <- lapply(igraph_list, function(x) blk_apply(x, concor_out, v))

  return(igraph_out)
}

concor_igraph_apply <- function(igraph_list, nsplit = 1) {
  adj_list <- lapply(igraph_list, function(x) igraph::get.adjacency(x, sparse = FALSE))

  concor_out <- suppressWarnings(concor(adj_list, p = nsplit))
  v <- paste("csplit", nsplit, sep = "")
  igraph_out <- lapply(igraph_list, function(x) blk_apply(x, concor_out, v))

  return(igraph_out)
}

concor_plot <- function(iobject, p = NULL) {
  split_name <- paste0("csplit", p)
  igraph::plot.igraph(iobject, vertex.color = igraph::vertex.attributes(iobject)[[split_name]],
                      vertex.label = NA, vertex.size = 5, edge.arrow.size = .3)
}

