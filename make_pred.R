library(ShortRead)
library(ape)
library(phangorn)
library(DECIPHER)
library(ggplot2)

nexp <- function(x) exp(1)^(-x)

get_node_label <- function(tree,idx){
  if (idx > length(tree$tip.label)){
    node <- tree$node.label[idx - length(tree$tip.label)]
  }else{
    node <- tree$tip.label[idx]
  }
  return(node)
}

predict_tip <- function(tree,node,mrra,mrra_idx,traits,node_dist){
  
  parent <- Ancestors(tree,node,type='parent')
  children <- Children(tree,parent)
  
  if (is.null(mrra)){
    
    ancestor_trait <- NULL
    ancestor_d <- NULL
    ancestor_w <- NULL
    
  }else{
    
    ancestor_trait <- traits[mrra,]
    ancestor_d <- node_dist[parent,mrra_idx]
    ancestor_w <- nexp(ancestor_d)
    
  }
  
  if (is.null(ancestor_trait)){
    
    pred <- NULL
    total_w <- NULL
    
  }else{
    
    pred <- ancestor_trait * ancestor_w
    total_w <- ancestor_w
    
  }
  
  for (i in children){
    
    child <- get_node_label(tree,i)
    
    if (child %in% rownames(traits)){
      
      child_trait <- traits[child,]
      child_parent_d <- node_dist[parent,i]
      child_parent_w <- nexp(child_parent_d)
      
      if (is.null(pred) && is.null(total_w)){
        
        total_w <- child_parent_w
        pred <- child_trait * child_parent_w
        
      }else{
        
        pred <- pred + child_trait * child_parent_w
        total_w <- total_w + child_parent_w
        
      }
      
    }
    
  }
  
  if (is.null(pred)) return(NA)
  
  pred <- pred/total_w
  
  return(pred)
  
}

create_alignemnt <- function(fasta, verbose=FALSE){
  
  if (verbose) cat('Performing alignment... ')
  seqs <- sread(fasta)
  names(seqs) <- as.character(id(fasta))
  
  alignment <- AlignSeqs(seqs,
                         processors=ncores,
                         useStructures=FALSE,
                         restrict=-1000,
                         verbose=verbose)
  
  alignment <- phyDat(as(alignment,'matrix'),type='DNA')
  if (verbose) cat('complete.\n')
  
  return(alignment)
  
}

create_tree <- function(alignment,fun_dist,fun_tree,ncores=1,verbose=FALSE){
  
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
  tree$node.label <- paste0('internal_node_',seq_along(unique(tree$edge[,1])))
  tree$edge.length[tree$edge.length < .0001] <- .0001
  
  attr(tree,'path') <- attr(fasta,'path')
  
  return(tree)
  
}

prune_tree_and_traits <- function(tree,img_x,otu_img){
  
  otu_img <- otu_img[!duplicated(otu_img$img_genome_id,fromLast=TRUE),]
  rownames(otu_img) <- otu_img$img_genome_id
  otu_img <- otu_img[otu_img$img_genome_id %in% img_x[,1],] #map must be first column
  otu_img <- otu_img[otu_img$X.gg_id %in% tree$tip.label,]
  
  tree_pruned <- drop.tip(tree,tip=tree$tip.label[!(tree$tip.label %in% otu_img$X.gg_id)])
  
  traits_pruned <- img_x[img_x[,1] %in% otu_img$img_genome_id,]
  rownames(traits_pruned) <- otu_img[as.character(traits_pruned[,1]),'X.gg_id']
  traits_pruned <- traits_pruned[,-1,drop=FALSE]
  
  return(list(tree=tree_pruned,traits=traits_pruned))
  
}


format_tree_and_traits <- function(tree,verbose=FALSE){
  
  if (verbose) cat('Loading tree and traits... ')
  img_ko <- read.delim(file.path(attr(tree,'path'),'IMG_ko_counts.tab'),sep='\t')
  img_16s <- read.delim(file.path(attr(tree,'path'),'IMG_16S_counts.tab'),sep='\t')
  otu_img <- read.delim(file.path(attr(tree,'path'),'gg_13_5_img.txt'),sep='\t')
  if (verbose) cat('complete.\n')
  
  
  if (verbose) cat('Formatting tree and traits... ')
  pruned_ko <- prune_tree_and_traits(tree,img_ko,otu_img)
  pruned_ko$tree$node.label[pruned_ko$tree$node.label == get_node_label(pruned_ko$tree,getRoot(pruned_ko$tree))] <- 'root'
  pruned_16s <- prune_tree_and_traits(tree,img_16s,otu_img)
  pruned_16s$tree$node.label[pruned_16s$tree$node.label == get_node_label(pruned_16s$tree,getRoot(pruned_16s$tree))] <- 'root'
  tree$node.label[tree$node.label == get_node_label(tree,getRoot(tree))] <- 'root'
  if (verbose) cat('complete.\n')
  
  out <- list(ko=pruned_ko,`16s`=pruned_16s,tree=tree)
  
  return(out)
  
}


asr <- function(tree,traits){
  
  traits <- traits[tree$tip.label,,drop=FALSE]
  
  reconstructions <- apply(traits,2,ace,tree,type='continuous',method='pic') 
  
  asr_traits <- lapply(seq_along(reconstructions),function(x) reconstructions[[x]]$ace)
  names(asr_traits) <- names(reconstructions)
  
  asr_traits <- do.call(cbind,asr_traits)
  rownames(asr_traits) <- tree$node.label
  
  return(asr_traits)
  
}

predict_traits <- function(traits,traits_asr,tree){
  
  traits <- traits[,colnames(traits_asr),drop=FALSE]
  traits <- as.matrix(rbind(traits_asr,traits))
  
  node_dist <- dist.nodes(tree)
  
  recons <- matrix(0.0,length(tree$tip.label),ncol(traits),dimnames=list(tree$tip.label,colnames(traits)))
  for (i in seq_along(tree$tip.label)){
    
    node_to_predict <- get_node_label(tree,i)
    
    cat(sprintf('%s: Tip %s... ',i,node_to_predict))
    
    ancestors <- Ancestors(tree,i,type='all')
    mrra <- NULL
    for (j in ancestors){
      ancestor <- get_node_label(tree,j)
      if (ancestor %in% rownames(traits)){
        mrra <- ancestor
        mrra_idx <- j
        break
      }
    }
    
    cat(sprintf('predicting with mrra %s\n',mrra))
    
    pred <- round(predict_tip(tree,i,mrra,mrra_idx,traits,node_dist))
    
    if (node_to_predict %in% rownames(traits)){
      recons[node_to_predict,] <- traits[node_to_predict,]
    }else{
      recons[node_to_predict,] <- pred
    }
    
  }
  
  return(recons)
  
}

create_ml_tree <- function(alignment,treeNJ,ncores=ncores,verbose=FALSE){
  
  if(verbose) cat('computes the likelihood of a phylogenetic tree... ')  
  tree_fit <- pml(treeNJ,data=alignment)
  if(verbose) cat('complete. \n')
  
  if(verbose) cat('Update the tree... ')
  fitGTR <- update(tree_fit,k=4,inv=0.2)
  if(verbose) cat('complete. \n')
  
  if(verbose) cat('computing the optim.pml ')
  fitGTR <- optim.pml(fitGTR,
                      model='F81',
                      optInv=TRUE,
                      optGamma=TRUE,
                      rearrangement='NNI',
                      control=pml.control(trace=1L))
  if(verbose) cat('complete \n')
  
  
  return(fitGTR$tree)
  
}



# read the alignment 
fasta <- readRDS('~/Documents/rcrust/alignment_rsv.rds')
#alignment <- create_alignemnt(fasta, verbose=TRUE)

#saveRDS(alignment,'~/Documents/rcrust/alignment_bench.rds')

alignment <- readRDS('~/Documents/rcrust/alignment_bench.rds')

tree <- create_tree(alignment,dist.ml,nj,ncores=50,verbose=TRUE)
# write.tree(tree,'~/Documents/rcrust/testing_workflow/rsv.tree')

####
ncores <- 6
treePML <- create_ml_tree(fasta,tree,ncores=ncores,verbose=TRUE)

tree_and_traits <- format_tree_and_traits(treePML,verbose=TRUE)

# get the ancestors
asr_traits_ko <- asr(tree_and_traits$ko$tree,tree_and_traits$ko$traits)
asr_traits_16s <- asr(tree_and_traits$`16s`$tree,tree_and_traits$`16s`$traits)

# save the asr and compare
saveRDS(asr_traits_ko,'~/Documents/rcrust/testing_workflow/asr_rsv_ko.rds')
saveRDS(asr_traits_16s,'~/Documents/rcrust/testing_workflow/asr_rsv_16s.rds')

reconstruction_ko <- predict_traits(tree_and_traits$ko$traits,asr_traits_ko,tree_and_traits$tree)
reconstruction_16s <- predict_traits(tree_and_traits$`16s`$traits,asr_traits_16s,tree_and_traits$tree)

saveRDS(reconstruction_ko,'~/Documents/rcrust/testing_workflow/reconstruction_ko.rds')
saveRDS(reconstruction_16s,'~/Documents/rcrust/testing_workflow/reconstruction_16s.rds')