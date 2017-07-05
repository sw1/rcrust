library(ShortRead)
library(ape)
library(phangorn)
library(DECIPHER)

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

create_fasta <- function(abundance,path,rows_are_taxa,verbose=FALSE){
  
  path_ref <- file.path(path,'gg_13_5.fasta.gz')
  path_img_ko <- file.path(path,'IMG_ko_counts.tab')
  path_otu_img <- file.path(path,'gg_13_5_img.txt')
  
  if (verbose) cat('Reading references... ')
  seqs_ref <- readFasta(path_ref)
  img_ko <- read.table(file=path_img_ko,sep='\t',row.names=1,header=TRUE,comment.char='')
  otu_img <- read.table(file=path_otu_img,sep='\t',header=TRUE,comment.char='')
  if (verbose) cat('complete.\n')
  
  if (verbose) cat('Preparing reference database... ')
  img_ko <- img_ko[rownames(img_ko) %in% otu_img$img_genome_id,]
  otu_img <- otu_img[otu_img$img_genome_id %in% rownames(img_ko),]
  
  seqs_ref <- seqs_ref[as.character(id(seqs_ref)) %in% otu_img$X.gg_id,]
  otus_ref <- as.character(id(seqs_ref))
  seqs_ref <- as.character(sread(seqs_ref))
  if (verbose)  cat('complete.\n')
  
  if (!rows_are_taxa) abundance <- t(abundance)
  
  if (verbose) cat('Creating fasta... ')
  seqs_sample <- rownames(abundance)
  rsvs_sample <- paste('rsv',1:nrow(abundance),sep='_')
  
  ids <- c(otus_ref,rsvs_sample)
  seqs <- c(seqs_ref,seqs_sample)
  fasta <- ShortRead(sread=DNAStringSet(seqs),id=BStringSet(ids))
  if (verbose) cat('complete.\n')
  
  attr(fasta,'path') <- path
  
  return(fasta)
  
}

create_tree <- function(fasta,fun_dist,fun_tree,ncores=1,verbose=FALSE){
  
  fun_dist <- match.fun(fun_dist)
  fun_tree <- match.fun(fun_tree)
  
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

abundance <- readRDS('~/rcrust/seqtab.rds')
path <- '~/rcrust/reference_files'

fasta <- create_fasta(abundance,path,FALSE,verbose=TRUE)
tree <- create_tree(fasta,dist.ml,nj,ncores=50,verbose=TRUE)
write.tree(tree,'~/rcrust/testing_workflow/rtree.tree')

tree_and_traits <- format_tree_and_traits(tree,verbose=TRUE)

asr_traits_ko <- asr(tree_and_traits$ko$tree,tree_and_traits$ko$traits)
asr_traits_16s <- asr(tree_and_traits$`16s`$tree,tree_and_traits$`16s`$traits)

reconstruction_ko <- predict_traits(tree_and_traits$ko$traits,asr_traits_ko,tree_and_traits$tree)
reconstruction_16s <- predict_traits(tree_and_traits$`16s`$traits,asr_traits_16s,tree_and_traits$tree)



# answer
final_picrust_table_ko <- read.delim('~/rcrust/testing_workflow/out/new_refs/ko_13_5_precalculated.tab.gz',sep='\t')
rownames(final_picrust_table_ko) <- final_picrust_table_ko$X.OTU_IDs
final_picrust_table_ko <- final_picrust_table_ko[,-1]

final_picrust_table_16s <- read.delim('~/rcrust/testing_workflow/out/new_refs/16S_13_5_precalculated.tab.gz',sep='\t')
rownames(final_picrust_table_16s) <- final_picrust_table_16s$X.OTU_IDs
final_picrust_table_16s <- final_picrust_table_16s[,-1,drop=FALSE]


# check
which(rowSums(reconstruction_ko) != rowSums(final_picrust_table_ko[rownames(reconstruction_ko),]))
which(reconstruction_16s != final_picrust_table_16s[rownames(reconstruction_16s),])


# ddna <- dist.dna(alignment_phydat)
# dml <- dist.ml(alignment_phydat)
# dmh <- dist.hamming(alignment_phydat)
# 
# tr1 <- NJ(dmh)
# tr2 <- nj(dmh)
# tr3 <- bionj(dmh)
# tr4 <- fastme.bal(dmh)
# tr5 <- fastme.ols(dmh)
# tr6 <- bionj(dmh)
# treedist(tr1,r2)

# ape::write.tree(tree,'~/rcrust/tree.tree')



# 
# writeFasta(fasta,file='~/rcrust/gg_13_5_study_db.fasta')
# write.table(as.character(id(fasta)),file='~/rcrust/traits_sample_filter.txt',sep='\t',quote=FALSE,col.names=FALSE,row.names=FALSE)
# 
# saveRDS(fasta,'~/rcrust/fasta.rds')


