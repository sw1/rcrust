library(ShortRead)
library(ape)
library(phangorn)
library(DECIPHER)

#args <- commandArgs(TRUE)
#print(args)

#if (length(args) != 3) stop('Provide 3 args: index, cores, seed')
#i <- as.integer(args[1])
#nc <- as.integer(args[2])
#seed <- as.integer(args[3])

#cat(sprintf('Running simulation for %s with %s cores, seed %s.\n',i,nc,seed))

nexp <- function(x) exp(1)^(-x)

get_node_label <- function(tree,idx){
  if (idx > length(tree$tip.label)){
    node <- tree$node.label[idx - length(tree$tip.label)]
  }else{
    node <- tree$tip.label[idx]
  }
  return(node)
}

create_tree_sim <- function(alignment,fun_dist,fun_tree,ncores=1,verbose=FALSE){
  #set.seed(seed)
    
  fun_dist <- match.fun(fun_dist) 
  fun_tree <- match.fun(fun_tree)
  

  if (verbose) cat('Computing pairwise distances... ')
  d <- fun_dist(alignment)
  if (verbose) cat('complete.\n')
  
  if (verbose) cat('Creating tree... ')
  tree <- fun_tree(d) 
  if (verbose) cat('complete.\n')
  
  if (verbose) cat('Formatting tree and traits... ')
  tree <- multi2di(tree,random=FALSE)
  # tree$node.label <- paste0('internal_node_',seq_along(unique(tree$edge[,1])))
  # tree$edge.length[tree$edge.length < .0001] <- .0001
  
  return(tree)
  
}

seed <- sample(seq_len(10000),1)
set.seed(seed)
ncores <- 6 #nc

seq_path <- '~/Documents/rcrust/seqtab.rds'
ref_path <- '~/Documents/rcrust/reference_files'
out_path <- '~/Documents/rcrust/testing_workflow/rtree_ml.tree'

# change dist.ml
# change nj

alignment <- readRDS(file.path(ref_path,'alignment.rds'))
tree <- create_tree_sim(alignment,dist.ml,nj,ncores=ncores,verbose=TRUE)
write.tree(tree,out_path)

# data(Laurasiatherian)
# tree <- create_tree_sim(Laurasiatherian,dist.hamming,nj,ncores=ncores,verbose=TRUE)
# cat(sprintf('Writing tree to %s.\n',out_path))
# write.tree(tree,out_path)
