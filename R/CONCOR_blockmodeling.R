
make_blk <- function(adj_list, splitn = 1) {
  #returns raw blockmodel data

  #run concor to get igraph splits, and shhh
  concor_out <- suppressWarnings(concor(adj_list, p = splitn))

  #check and reorder the concor output to be the same as the vertex list
  concor_order <- match(colnames(adj_list[[1]]), concor_out$vertex)
  block_ordered <- concor_out$block[concor_order]

  # ##prepare to attach sna (needed for next part) so it can later be detached easily
  # #record the packages that existed at the start
  # pre <- names(sessionInfo()$otherPkgs)
  # #get rid of all nondefault packages
  # if (length(pre) != 0) {
  #   invisible(lapply(paste('package:', names(sessionInfo()$otherPkgs), sep = ""), detach,
  #                    character.only = TRUE, unload = TRUE))
  # }
  # #add that asshole of a package nice and quietly
  # suppressMessages(library(sna))

  #make blokmodels (aka the thing all this has lead to)
  blockmodel_list <- lapply(adj_list, function(x) sna::blockmodel(as.matrix(x), block_ordered))

  # ##detach sna and reattach whatever else
  # suppressMessages(invisible(lapply(paste('package:', names(sessionInfo()$otherPkgs), sep = ""),
  #                                   detach, character.only = TRUE, unload = TRUE)))
  # suppressMessages(invisible(lapply(pre, library, character.only = TRUE)))

  #return raw blockmodel data
  return(blockmodel_list)
}

edge_dens <- function(adj_mat) {
  adj_mat[adj_mat > 0] <- 1
  a <- sum(adj_mat)
  m <- length(adj_mat) - sqrt(length(adj_mat))
  d <- a / m
  return(d)
}

make_reduced <- function(adj_list, splitn = 1, weighted = FALSE) {
  #returns list of reduced matrixes one for each relation in $red_mat
  #also returns the cutoff densities in $dens
  #inputs are list of adj matrixes what split you want and weather weighted or not (True is weighted)

  #get the raw blockmodel results
  blk_out = make_blk(adj_list, splitn)

  #check density of the adjacency matrix
  dens_vec <- sapply(adj_list, function(x) edge_dens(x))

  #pull out the density matrix of blk_out
  d <- lapply(blk_out, function(x) x[[5]])

  #play with the matrix for later use
  mat_return <- vector("list", length = length(dens_vec))

  for (i in 1:length(dens_vec)) {
    temp1 <- d[[i]]
    #check for and make NaN to 0
    temp1[is.nan(temp1)] <- 0
    #set blocks with density less than overall density to 0
    temp1[temp1 < dens_vec[[i]]] <- 0
    if (!weighted) {
      temp1[temp1 > 0] <- 1
      mat_return[[i]] <- temp1
    }
    if (weighted) {
      #find minimum density and scale all by that
      min <- min(temp1[temp1 > 0])
      mat_return[[i]] <- as.matrix(temp1 / min)
      #scale all densities to 10 or less
      while (max(mat_return[[i]]) > 20) {
        mat_return[[i]] <- mat_return[[i]] / 1.05
      }
      while (max(mat_return[[i]]) < 18) {
        mat_return[[i]] <- mat_return[[i]] * 1.05
      }
    }
  }
  return_list <- list()
  return_list$reduced_mat <- mat_return
  return_list$dens <- dens_vec
  return(return_list)
}

plot_blk <- function (x, ...) {
  #slightly edited version of the function from SNA, plots as square and slightly changed labeling also attached SNA and later removed
  #why is it convention to not have comments

  # ##attach SNA for later use and record packages that existed beforehand
  # #record the packages that existed at the start
  # pre <- names(sessionInfo()$otherPkgs)
  # #get rid of all nondefault packages
  # invisible(lapply(paste('package:', names(sessionInfo()$otherPkgs), sep = ""), detach,
  #                  character.only = TRUE, unload = TRUE))
  # #add sna
  # suppressMessages(library(sna))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  n <- dim(x$blocked.data)[2]
  m <- sna::stackcount(x$blocked.data)
  if (!is.null(x$plabels))
    plab <- x$plabels
  else plab <- (1:n)[x$order.vector]
  glab <- ""
  par(mfrow = c(floor(sqrt(m)), ceiling(m/floor(sqrt(m)))))
  if (m > 1)
    for (i in 1:m) {
      sna::plot.sociomatrix(x$blocked.data[i, , ], labels = list(plab, plab),
                       main = glab[i], drawlines = FALSE, asp = 1)

      for (j in 2:n) if (x$block.membership[j] != x$block.membership[j - 1])
        abline(v = j - 0.5, h = j - 0.5, lty = 3)
    }
  else {
    sna::plot.sociomatrix(x$blocked.data, labels = list(plab, plab),
                     main = glab[1], drawlines = FALSE, asp = 1)

    for (j in 2:n) if (x$block.membership[j] != x$block.membership[j - 1])
      abline(v = j - 0.5, h = j - 0.5, lty = 3)
  }

  # ##detach sna and reattach whatever else
  # invisible(lapply(paste('package:', names(sessionInfo()$otherPkgs), sep = ""),
  #                  detach, character.only = TRUE, unload = TRUE))
  # suppressMessages(invisible(lapply(pre, library, character.only = TRUE)))
}

plot_blk_labeless <- function(bm) {
  #get rid of those pesky labels
  bm$plabels <- rep("", length(bm$plabels))
  bm$glabels <- ""
  plot_blk(bm)
}


#Igraph needed beyond this point
make_reduced_igraph <- function(reduced_mat) {
  w <- NULL
  if (any(reduced_mat > 1)) {
    w <- TRUE
  }
  iplotty <- igraph::graph_from_adjacency_matrix(reduced_mat, mode = "directed", weighted = w)
  return(iplotty)
}

plot_red_weighted <- function(blk) {
  #plots just weighted blockmodel network
  #plots colors based on order of blocks (I think this should match those on the networks)
  igraph::plot.igraph(blk, vertex.color = c(1:length(igraph::vertex.attributes(blk)[[1]])), vertex.label = NA,
       edge.width = (E(blk)$weight/3), edge.arrow.size = (igraph::E(blk)$weight/15), vertex.size = 25)
}

plot_red_unweighted <- function(blk) {
  #plots just weighted blockmodel network
  #plots colors based on order of blocks (I think this should match those on the networks)
  igraph::plot.igraph(blk, vertex.color = c(1:length(igraph::vertex.attributes(blk)[[1]])), vertex.label = NA,
       edge.arrow.size = .6, vertex.size = 25)
}

