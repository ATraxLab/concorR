#' @export
make_blk <- function(adj_list, nsplit = 1) {
  concor_out <- suppressWarnings(concor(adj_list, nsplit))

  concor_order <- match(colnames(adj_list[[1]]), concor_out$vertex)
  block_ordered <- concor_out$block[concor_order]

  blockmodel_list <- lapply(adj_list,
                            function(x) sna::blockmodel(as.matrix(x),
                                                        block_ordered))

  return(blockmodel_list)
}

.edge_dens <- function(adj_mat) {
  adj_mat[adj_mat > 0] <- 1

  a <- sum(adj_mat)
  m <- length(adj_mat) - sqrt(length(adj_mat))
  d <- a / m
  return(d)
}

.scaledDegree <-  function(adj_mat){
    adj_mat[adj_mat > 0] <- 1

    avgOutDegree = sum(adj_mat)/nrow(adj_mat)
    maxOutDegree = max(rowSums(adj_mat))
    scaledDegree = ifelse(maxOutDegree > 0,avgOutDegree/maxOutDegree, 0)
    return(scaledDegree)
}

#' @export
make_reduced <- function(adj_list, nsplit = 1, stat='density') {
  if(stat=='density'){  
    blk_out = make_blk(adj_list, nsplit)
    dens_vec <- sapply(adj_list, function(x) .edge_dens(x))
    d <- lapply(blk_out, function(x) x[[5]])
    mat_return <- vector("list", length = length(dens_vec))

    for (i in 1:length(dens_vec)) {
      temp1 <- d[[i]]
      temp1[is.nan(temp1)] <- 0
      temp1[temp1 < dens_vec[[i]]] <- 0
      temp1[temp1 > 0] <- 1
      mat_return[[i]] <- temp1
    }
  
    return_list <- list()
    return_list$reduced_mat <- mat_return
    return_list$dens <- dens_vec
    return(return_list)
  }else if(stat=='degree'){
    blk_out = make_blk(adj_list, nsplit)
    outdegree = lapply(adj_list, function(x) .scaledDegree(x))
    mat_return <- vector("list", length = length(outdegree))
    
    for(i in 1:length(outdegree)){ 
      this_adj_mat = adj_list[[i]]
      thisBlk = blk_out[[i]]
      members = thisBlk$block.membership[sort(thisBlk$order.vector,index.return=TRUE)$ix]
      nb = max(members)
      reduced_degree = matrix(0, nrow = nb, ncol = nb)
      rownames(reduced_degree) = paste("Block",1:nb)
      colnames(reduced_degree) = paste("Block",1:nb)
      for(j in 1:nb){
        nRows = sum(j==members)
        for(k in 1:nb){
          nCols = sum(k==members)
          if(nRows==1){
            if(nCols==1){
              blk_adj_mat = this_adj_mat[j==members, k==members] 
              outDeg = ifelse(blk_adj_mat>0,1,0) 
            }else{
              blk_adj_mat = this_adj_mat[j==members, k==members] 
              blk_adj_mat = matrix(blk_adj_mat,nrow=1)
              outDeg = .scaledDegree(blk_adj_mat)
            }
          }else{
            if(nCols==1){
              blk_adj_mat = this_adj_mat[j==members, k==members]
              blk_adj_mat = matrix(blk_adj_mat,ncol=1)
            }else{
              blk_adj_mat = this_adj_mat[j==members, k==members]
            }
            outDeg = .scaledDegree(blk_adj_mat)
          }
          reduced_degree[j,k] = outDeg
        }
      }
      temp1 <- reduced_degree
      temp1[is.nan(temp1)] <- 0
      temp1[temp1 < outdegree[[i]]] <- 0
      temp1[temp1 > 0] <- 1
      mat_return[[i]] <- temp1
    }
                                                                
    return_list <- list()
    return_list$reduced_mat <- mat_return
    return_list$deg <- outdegree
    return(return_list)
  }else{
    stop('Statistics implemented for determine edges in reduced networks are only 
         density and degree.')
  }
}

#' @export
plot_blk <- function (x, labels = FALSE, ...) {
  # Adapted from sna::plot.blockmodel():
  # Carter T. Butts (2019). sna: Tools for Social Network Analysis.
  # R package version 2.5, licensed under GPL (>= 2).
  # https://CRAN.R-project.org/package=sna
  if (!labels) {
    x$plabels <- rep("", length(x$plabels))
    x$glabels <- ""
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  n <- dim(x$blocked.data)[2]
  m <- sna::stackcount(x$blocked.data)
  if (!is.null(x$plabels))
    plab <- x$plabels
  else plab <- (1:n)[x$order.vector]
  glab <- ""
  graphics::par(mfrow = c(floor(sqrt(m)), ceiling(m/floor(sqrt(m)))))
  if (m > 1)
    for (i in 1:m) {
      sna::plot.sociomatrix(x$blocked.data[i, , ], labels = list(plab, plab),
                       main = glab[i], drawlines = FALSE, asp = 1)

      for (j in 2:n) if (x$block.membership[j] != x$block.membership[j - 1])
        graphics::abline(v = j - 0.5, h = j - 0.5, lty = 3)
    }
  else {
    sna::plot.sociomatrix(x$blocked.data, labels = list(plab, plab),
                     main = glab[1], drawlines = FALSE, asp = 1)

    for (j in 2:n) if (x$block.membership[j] != x$block.membership[j - 1])
      graphics::abline(v = j - 0.5, h = j - 0.5, lty = 3)
  }
}

#' @export
make_reduced_igraph <- function(reduced_mat) {
  reduced_graph <- igraph::graph_from_adjacency_matrix(reduced_mat,
                                                       mode = "directed")
  return(reduced_graph)
}


#' @export
plot_reduced <- function(iobject,main='') {
    vcolors  <- viridis::viridis(length(igraph::vertex_attr(iobject)$name))
    igraph::plot.igraph(iobject, vertex.color = vcolors, vertex.label = NA,
                      edge.arrow.size = .6, vertex.size = 25,main = main)
}

.block_edge_dens  <- function(adj_mat){
  adj_mat[adj_mat > 0]  <- 1
  
  a <- sum(adj_mat)
  m <- length(adj_mat)
  d <- a / m
  return(d)
}    

#' @export
make_reduced_from_partition <- function(adj_mat, partition, stat='density') {
  if(stat=='density'){  
    dens <- .edge_dens(adj_mat)
    
    nb = max(partition)
    reduced_den = matrix(0, nrow = nb, ncol = nb)
    rownames(reduced_den) = paste("Block",1:nb)
    colnames(reduced_den) = paste("Block",1:nb)
    for(j in 1:nb){
      nRows = sum(j==partition)
      for(k in 1:nb){
        nCols = sum(k==partition)
        if(nRows==1){
          if(nCols==1){
            blk_adj_mat = adj_mat[j==partition, k==partition] 
            d = ifelse(blk_adj_mat>0,1,0) 
          }else{
            blk_adj_mat = adj_mat[j==partition, k==partition] 
            blk_adj_mat = matrix(blk_adj_mat,nrow=1)
            d = .block_edge_dens(blk_adj_mat)
          }
        }else{
          if(nCols==1){
            blk_adj_mat = adj_mat[j==partition, k==partition]
            blk_adj_mat = matrix(blk_adj_mat,ncol=1)
          }else{
            blk_adj_mat = adj_mat[j==partition, k==partition]
          }
          d = ifelse(j==k,.edge_dens(blk_adj_mat),
                     .block_edge_dens(blk_adj_mat))
        }
        reduced_den[j,k] = d
      }
    }
    reduced_den[is.nan(reduced_den)] <- 0
    reduced_den[reduced_den < dens] <- 0
    reduced_den[reduced_den > 0] <- 1
    
    return_list <- list()
    return_list$reduced_mat <- reduced_den
    return_list$dens <- dens
    return(return_list)
  }else if(stat=='degree'){
    outdegree = .scaledDegree(adj_mat)
    
    nb = max(partition)
    reduced_degree = matrix(0, nrow = nb, ncol = nb)
    rownames(reduced_degree) = paste("Block",1:nb)
    colnames(reduced_degree) = paste("Block",1:nb)
    for(j in 1:nb){
      nRows = sum(j==partition)
      for(k in 1:nb){
        nCols = sum(k==partition)
        if(nRows==1){
          if(nCols==1){
            blk_adj_mat = adj_mat[j==partition, k==partition] 
            outDeg = ifelse(blk_adj_mat>0,1,0) 
          }else{
            blk_adj_mat = adj_mat[j==partition, k==partition] 
            blk_adj_mat = matrix(blk_adj_mat,nrow=1)
            outDeg = .scaledDegree(blk_adj_mat)
          }
        }else{
          if(nCols==1){
            blk_adj_mat = adj_mat[j==partition, k==partition]
            blk_adj_mat = matrix(blk_adj_mat,ncol=1)
          }else{
            blk_adj_mat = adj_mat[j==partition, k==partition]
          }
          outDeg = .scaledDegree(blk_adj_mat)
        }
        reduced_degree[j,k] = outDeg
      }
    }
    reduced_degree[is.nan(reduced_degree)] <- 0
    reduced_degree[reduced_degree < outdegree] <- 0
    reduced_degree[reduced_degree > 0] <- 1
    
    return_list <- list()
    return_list$reduced_mat <- reduced_degree
    return_list$deg <- outdegree
    return(return_list)
  }else{
    stop('Statistics implemented for determining edges in reduced networks are only 
         density and degree.')
  }
}

