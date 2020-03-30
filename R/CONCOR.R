#CONCOR and needed functions
#created 11/26/18
#By Tyme Suda
#updated 8/20/18

#All needed functions to run the concor function for an arbituary number of
#splits (equal or less than the maximimum for said data)
#Works for single or maultiple relations formated as square matrixies in a simmple list
#*should* work for weighted networks would be nice to have known resaults to compare to


concor1 <- function(m_stack, cutoff = .9999999, max_iter = 50) {
  #does concor once
  #needs the matrix stack as an input
  #outputs a matrix
  if (ncol(m_stack) < 2)
    stop("Too few columns to partition.")
  cor_i <- cor(m_stack, use  =  "pairwise.complete.obs")
  iter<-0
  while (any(abs(cor_i) <= cutoff, na.rm = TRUE) & iter <= max_iter) {
    cor_i <- cor(cor_i, use = "pairwise.complete.obs")
    iter <- iter + 1
  }
  cor_i[cor_i < 0] <- -1
  cor_i[cor_i > 0] <- 1
  return(cor_i)
}

.val_diag <- function(m, value = NA) {
  #make  the diaganlas of a matrix be value
  #value defaults to NA
  #return matrix
  diag(m) <- value
  return(m)
}

.isolates_col <- function(m) {
  #taake in matrix return return col names of isolate columns
  #diaganals must be stet to zero beforhand
  #check to see what nodes are isolates
  #return vector of isolate names
  bob <- raw()
  for (i in 1:length(colnames(m))) {
    if (all(m[, i] == 0)) {
      bob <- c(bob, colnames(m)[i])
    }
  }
  return(bob)
}

.stack_mat <- function(m_list) {
  #makes stack of matrixes followed by their transposes for calling in concor
  #input a list of matrixes, output a tall matrix
  mt_list <- lapply(m_list, t)
  m_mt_list <- c(m_list, mt_list)
  mat_stack <- do.call("rbind", m_mt_list)
}

.make_order <- function(order_list) {
  #conbines suborders into an overall order for applying to the previously orderd matrix inputs
  #takes in a list of orders, returns a order vector
  if (!is.list(order_list)) {
    stop("not a list")
  }
  num_list <- length(order_list)
  l_sub <- 0
  for (i in 1:num_list) {
    order_list[[i]] <- order_list[[i]] + n
    n <- n + length(order_list[[i]])
  }
  order <- unlist(order_list)
  return(order)
}

.order_apply <- function(order, mat) {
  #orders rows and collumns of mat for use via lapply
  #input is the order being applied followed by the matrix it is applied to
  #output is orderd matrix
  m1 <- mat[order, order]
  return(m1)
}

.make_sub_boolean <- function(cor_matrixies_orderd) {
  #make boolean of first row of cor_mat_orderd where -1=>FALSE 1=>TRUE
  #for use with lapply, input matrix, output boolean vector
  group <- cor_matrixies_orderd[, 1] > 0
  return(group)
}

.make_big_booleans <- function(bool_list) {
  #makes 2 booleans (false=false, and false=true) that can be applied to the total matrix from each input boolean in bool_list
  #input is list of booleans
  #output is list double length of input list of booleans as long as all logical expersions in input
  if (!is.list(bool_list)) {
    stop("not a list")
  }
  bool_num <- length(bool_list)
  tot_length <- length(unlist(bool_list))
  #setup list of 2*boo_length false logical vectors each as long as tot_lenght
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
  #apply a boolean to a matrix, removing false collumns, check boolean correct size
  #input boolean vector, matrix
  #output matrix
  if (ncol(mat_stack) != length(boolean)) {
    stop("boolean of wrong size")
  }
  stack_shrunck <- mat_stack[, boolean, drop = FALSE]
  return(stack_shrunck)
}

.block_names <- function(mat_list) {
  lapply(seq_along(mat_list), function(x) data.frame(block = x,
                                                     vertex = colnames(mat_list[[x]]),
                                                     stringsAsFactors = FALSE))
}

concor <- function(m0, cutoff = .9999999, max_iter = 50, p = 1) {
  #Inputed m0 must be a list of matrixes WITH COL/ROW NAMES
  #outpust vectors of collumn names of final grouping

  #initializing varibles
  #set diaganals of each relation to zero for isolate surch
  mi <- lapply(m0, function(x) .val_diag(x, 0))
  miso <- mi

  if (length(.isolates_col(.stack_mat(miso))) > 0) {
    #remove isolates
    #make a stack of matrix and transposes to check overall isolates over
    iso_stack <- .stack_mat(miso)
    #check number of relations
    num_relat <- length(miso)
    #get overall isolates of all relations
    isolist <- .isolates_col(iso_stack)
    #make vector of names to compare to and make boolean
    mi_names <- colnames(miso[[1]])
    #make boolean of which are isolates
    #true is isloate, flase not isolate
    iso_bool <- mi_names %in% isolist
    #make list of matrixes without the isolates (apply iso_bool to each relation)
    for (i in 1:num_relat) {
      #make matrix with no isolates
      mi[[i]] <- miso[[i]][!iso_bool, !iso_bool, drop = FALSE]
    }
    #make matrix of just isolates
    m_iso <- m0[[1]][iso_bool, iso_bool, drop = FALSE]
  }


  #set diaganals to NA
  mi <- lapply(mi, function(x) .val_diag(x, NA))

  stack_list <- list(.stack_mat(mi))

  for (i in 1:p) {
    #apply concor1 to each matrix stack (should have 2^(i-1) mat stacks/elements of output list)
    concored <- lapply(stack_list, function(x) concor1(x))
    #make orders sub orders
    order_list <- lapply(concored, function(x) order(x[, 1]))
    #apply each suborder corisponding matrix in concored
    for (j in 1:(2^(i-1))) {
      concored[[j]] <- .order_apply(order_list[[j]], concored[[j]])
    }
    #conbine suborders
    order_list <- .make_order(order_list)
    #get sub booleans for next run
    bool_list <- lapply(concored, function(x) .make_sub_boolean(x))
    #make big booleans that can be applied
    bool_list <- .make_big_booleans(bool_list)
    #apply order to input matrixes
    mi <- lapply(mi, function(x) .order_apply(order_list, x))
    #make fat matrix stack
    stack_list <- .stack_mat(mi)
    #make list of skinny matrix stacks for next itteration (apply booleans)
    stack_list <- lapply(bool_list, function(x) .boolean_apply(x, stack_list))
  }
  #tack on the isolate group to stack list to run blocknames on
  mats_groups <- stack_list
  if (exists("m_iso")) {
    mats_groups[[length(stack_list)+1]] <- m_iso
  }
  #run blocknames to produce/format groups
  groups <- do.call(rbind, .block_names(mats_groups))
  groups[match(rownames(m0[[1]]), groups$vertex), ]

  return(groups)
}

