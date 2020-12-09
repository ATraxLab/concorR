#' @export
#' @importFrom stats cor
concor1 <- function(m_stack, cutoff = .9999999, max_iter = 50) {
  if (ncol(m_stack) < 2){
    return(cor(m_stack, use  =  "pairwise.complete.obs"))
  }

  cor_i <- cor(m_stack, use  =  "pairwise.complete.obs")
  iter <- 0

  while (any(abs(cor_i) <= cutoff, na.rm = TRUE) & iter <= max_iter) {
    cor_i <- cor(cor_i, use = "pairwise.complete.obs")
    iter <- iter + 1
  }
  cor_i[cor_i < 0] <- -1
  cor_i[cor_i > 0] <- 1

  return(cor_i)
}

.name <- function(mat) {
  a <- 1:nrow(mat)
  vnames <- sprintf("v%03d", a)
  colnames(mat) <- vnames
  rownames(mat) <- vnames
  return(mat)
}

.igname <- function(ig) {
  a <- 1:igraph::vcount(ig)
  vnames <- sprintf("v%03d", a)
  ig <- igraph::set_vertex_attr(ig, "name", value = vnames)
  return(ig)
}

.concor_validitycheck <- function(m_list) {
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

.concor_igraph_validitycheck <- function(ig_list) {
  b <- sapply(ig_list, function(x) !(igraph::is_named(x)))
  if (all(b)) {
    warning("node names don't exist\nAdding default node names\n")
    ig_list <- lapply(ig_list, function(x) .igname(x))
    b <- sapply(ig_list, function(x) !(igraph::is_named(x)))
  }
  if (any(b)) {
    stop("Node name mismatch")
  }

  return(ig_list)
}

.val_diag <- function(m, value = NA) {
  diag(m) <- value
  return(m)
}

.isolates_col <- function(m) {
  isolates <- raw()
  for (i in 1:length(colnames(m))) {
    if (all(m[, i] == 0)) {
      isolates <- c(isolates, colnames(m)[i])
    }
  }
  return(isolates)
}

.stack_mat <- function(m_list) {
  mt_list <- lapply(m_list, t)
  m_mt_list <- c(m_list, mt_list)
  mat_stack <- do.call("rbind", m_mt_list)
}

.make_order <- function(order_list) {
  if (!is.list(order_list)) {
    stop("not a list")
  }
  num_list <- length(order_list)
  n <- 0
  for (i in 1:num_list) {
    order_list[[i]] <- order_list[[i]] + n
    n <- n + length(order_list[[i]])
  }
  order <- unlist(order_list)
  return(order)
}

.order_apply <- function(order, mat) {
  if (length(order) == 1) {
    return(mat)
  }
  m1 <- mat[order, order]
  return(m1)
}

.make_sub_boolean <- function(cor_matrixies_ordered) {
  group <- cor_matrixies_ordered[, 1] > 0
  return(group)
}

.make_big_booleans <- function(bool_list) {
  if (!is.list(bool_list)) {
    stop("not a list")
  }
  bool_num <- length(bool_list)
  tot_length <- length(unlist(bool_list))
  booleans_out <- rep(list(vector("logical", tot_length)), 2 * bool_num)
  a <- 1
  for (i in 1:bool_num) {
    for(j in 1:length(bool_list[[i]])) {
      booleans_out[[2*i-1]][a] <- bool_list[[i]][j]
      booleans_out[[2*i]][a] <- !bool_list[[i]][j]
      a <- a + 1
    }
  }
  return(booleans_out)
}

.boolean_apply <- function(boolean, mat_stack) {
  if (ncol(mat_stack) != length(boolean)) {
    stop("boolean of wrong size")
  }
  stack_shrunck <- mat_stack[, boolean, drop = FALSE]
  return(stack_shrunck)
}

.block_names <- function(mat_list) {
  lapply(seq_along(mat_list),
         function(x) data.frame(block = x, vertex = colnames(mat_list[[x]]),
                                stringsAsFactors = FALSE))
}

#' @export
concor <- function(m_list, nsplit = 1, self_ties = FALSE, cutoff = .9999999,
                   max_iter = 50) {
  m_list <- .concor_validitycheck(m_list)
  mi <- m_list
  if (all(sapply(mi, function(x) all(is.na(diag(x)))))) {
    mi <- lapply(mi, function(x) .val_diag(x, 0))
  }
  miso <- mi

  if (length(.isolates_col(.stack_mat(miso))) > 0) {
    iso_stack <- .stack_mat(miso)
    num_relat <- length(miso)
    isolist <- .isolates_col(iso_stack)
    mi_names <- colnames(miso[[1]])
    iso_bool <- mi_names %in% isolist
    for (i in 1:num_relat) {
      mi[[i]] <- miso[[i]][!iso_bool, !iso_bool, drop = FALSE]
    }
    m_iso <- m_list[[1]][iso_bool, iso_bool, drop = FALSE]
  }

  if (!self_ties) {
    mi <- lapply(mi, function(x) .val_diag(x, NA))
  }
  stack_list <- list(.stack_mat(mi))
  stop_check <- list()

  for (i in 1:nsplit) {
    concored <- lapply(stack_list, function(x) concor1(x))
    order_list <- lapply(concored, function(x) order(x[, 1]))
    for (j in 1:length(order_list)) {
      concored[[j]] <- .order_apply(order_list[[j]], concored[[j]])
    }
    order_list <- .make_order(order_list)
    bool_list <- lapply(concored, function(x) .make_sub_boolean(x))
    bool_list <- .make_big_booleans(bool_list)
    mi <- lapply(mi, function(x) .order_apply(order_list, x))
    stack_list <- .stack_mat(mi)
    stack_list <- lapply(bool_list, function(x) .boolean_apply(x, stack_list))

    is_empty <- sapply(stack_list, function(x) ncol(x) != 0)
    stack_list <- stack_list[is_empty]

    if (identical(stop_check, stack_list)) {
      warning(paste("split", nsplit,
                    "was the same as split",  i - 1, "\n stopping"))
      break
    }

    stop_check <- stack_list
  }

  mats_groups <- stack_list
  if (exists("m_iso")) {
    mats_groups[[length( mats_groups) + 1]] <- m_iso
  }
  groups <- do.call(rbind, .block_names(mats_groups))
  groups[match(rownames(m_list[[1]]), groups$vertex), ]

  return(groups)
}

