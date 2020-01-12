#CONCOR and needed functions
#created 11/26/18
#By Tyme Suda
#updated 8/20/18

#All needed functions to run the concor function for an arbituary number of
#splits (equal or less than the maximimum for said data)
#Works for single or maultiple relations formated as square matrixies in a simmple list
#*should* work for weighted networks would be nice to have known resaults to compare to


library(igraph)

concor1=function(m_stack, cutoff=.9999999, max_iter=50){
  #does concor once
  #needs the matrix stack as an input
  #outputs a matrix
  if (ncol(m_stack) < 2)
    stop("Too few columns to partition.")
  cor_i=cor(m_stack, use = "pairwise.complete.obs")
  iter=0
  while (any(abs(cor_i) <= cutoff, na.rm=TRUE) & iter <= max_iter){
    cor_i=cor(cor_i, use = "pairwise.complete.obs")
    iter=iter+1
  }
  cor_i[cor_i<0]=-1
  cor_i[cor_i>0]=1
  return(cor_i)
}

val.diag =function(m, value=NA){
  #make  the diaganlas of a matrix be value
  #value defaults to NA
  #return matrix
  diag(m) <- value
  return(m)
}

isolates.col = function(m){
  #taake in matrix return return col names of isolate columns
  #diaganals must be stet to zero beforhand
  #check to see what nodes are isolates
  #return vector of isolate names
  bob=raw()
  for (i in 1:length(colnames(m))) {
    if (all(m[,i]==0)) {
      bob=c(bob,colnames(m)[i])
    }
  }
  return(bob)
}

stack_mat=function(m_list){
  #makes stack of matrixes followed by their transposes for calling in concor
  #input a list of matrixes, output a tall matrix
  mt_list=lapply(m_list, t)
  m_mt_list=c(m_list,mt_list)
  mat_stack=do.call("rbind", m_mt_list)
}

make_order=function(order_list){
  #conbines suborders into an overall order for applying to the previously orderd matrix inputs
  #takes in a list of orders, returns a order vector
  if (!is.list(order_list)){
    stop("not a list")
  }
  l_length=length(order_list)
  l_sub=0
  for (i in 1:l_length){
    order_list[[i]]=order_list[[i]]+l_sub
    l_sub=l_sub+length(order_list[[i]])
  }
  order=unlist(order_list)
  return(order)
}

order_apply=function(order,mat){
  #orders rows and collumns of mat for use via lapply
  #input is the order being applied followed by the matrix it is applied to
  #output is orderd matrix
  m1=mat[order,order]
  return(m1)
}

make_sub_boolian=function(cor_mat_orderd){
  #make boolian of first row of cor_mat_orderd where -1=>FALSE 1=>TRUE
  #for use with lapply, input matrix, output boolian vector
  group=cor_mat_orderd[, 1] > 0
  return(group)
}

make_big_boolians=function(bool_list){
  #makes 2 boolians (false=false, and false=true) that can be applied to the total matrix from each input boolian in bool_list
  #input is list of boolians
  #output is list double length of input list of boolians as long as all logical expersions in input
  if (!is.list(bool_list)){
    stop("not a list")
  }
  boo_num=length(bool_list)
  tot_length=length(unlist(bool_list))
  #setup list of 2*boo_length false logical vectors each as long as tot_lenght
  boolians_out=rep(list(vector("logical", tot_length)),2*boo_num)
  a=1
  for (i in 1:boo_num){
    for(j in 1:length(bool_list[[i]])){
      boolians_out[[2*i-1]][a]=bool_list[[i]][j]
      boolians_out[[2*i]][a]=!bool_list[[i]][j]
      a=a+1
    }
  }
  return(boolians_out)
}

boolian_apply=function(boolian, mat_stack){
  #apply a boolian to a matrix, removing false collumns, check boolian correct size
  #input boolian vector, matrix
  #output matrix
  if (ncol(mat_stack) != length(boolian)){
    stop("Boolian wrong size")
  }
  stack_shrunck=mat_stack[,boolian, drop=FALSE]
  return(stack_shrunck)
}

block_names=function(mat_list){
  lapply(seq_along(mat_list), function(x) data.frame(block = x,
                                                     vertex = colnames(mat_list[[x]]), stringsAsFactors = FALSE))
}

concor=function(m0, cutoff=.9999999, max_iter=50, p=1){
  #Inputed m0 must be a list of matrixes WITH COL/ROW NAMES
  #outpust vectors of collumn names of final grouping

  #initializing varibles
  #set diaganals of each relation to zero for isolate surch
  mi=lapply(m0, function(x) val.diag(x,0))
  miso=mi

  if (length(isolates.col(stack_mat(miso)))>0) {
    #remove isolates
    #make a stack of matrix and transposes to check overall isolates over
    iso.stack=stack_mat(miso)
    #check number of relations
    num.relat=length(miso)
    #get overall isolates of all relations
    isolist=isolates.col(iso.stack)
    #make vector of names to compare to and make boolian
    mi.names=colnames(miso[[1]])
    #make boolian of which are isolates
    #true is isloate, flase not isolate
    iso.bool = mi.names %in% isolist
    #make list of matrixes without the isolates (apply iso.bool to each relation)
    for (i in 1:num.relat) {
      #make matrix with no isolates
      mi[[i]]=miso[[i]][!iso.bool,!iso.bool, drop=FALSE]
    }
    #make matrix of just isolates
    m.iso=m0[[1]][iso.bool,iso.bool, drop=FALSE]
  }


  #set diaganals to NA
  mi=lapply(mi, function(x) val.diag(x, NA))

  stack_list=list(stack_mat(mi))

  for (i in 1:p){
    #apply concor1 to each matrix stack (should have 2^(i-1) mat stacks/elements of output list)
    concored=lapply(stack_list, function(x) concor1(x, cutoff=.9999999, max_iter=50))
    #make orders sub orders
    order_list=lapply(concored, function(x) order(x[,1]))
    #apply each suborder corisponding matrix in concored
    for (j in 1:(2^(i-1))){
      concored[[j]]=order_apply(order_list[[j]], concored[[j]])
    }
    #conbine suborders
    order_list=make_order(order_list)
    #get sub boolians for next run
    bool_list=lapply(concored, function(x) make_sub_boolian(x))
    #make big boolians that can be applied
    bool_list=make_big_boolians(bool_list)
    #apply order to input matrixes
    mi=lapply(mi, function(x) order_apply(order_list, x))
    #make fat matrix stack
    stack_list=stack_mat(mi)
    #make list of skinny matrix stacks for next itteration (apply boolians)
    stack_list=lapply(bool_list, function(x) boolian_apply(x, stack_list))
  }
  #tack on the isolate group to stack list to run blocknames on
  mats.groups= stack_list
  if (exists("m.iso")) {
    mats.groups[[length(stack_list)+1]]=m.iso
  }
  #run blocknames to produce/format groups
  groups <- do.call(rbind, block_names(mats.groups))
  groups[match(rownames(m0[[1]]), groups$vertex), ]

  return(groups)
}

