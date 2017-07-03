library(ShortRead)
library(ape)
library(phangorn)
library(DECIPHER)

nexp <- function(x) exp(1)^(-x) # -1x

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
    
    # child <- tree$tip.label[i]
    child <- get_node_label(tree,i)
    
    if (child %in% rownames(traits)){
      
      child_trait <- traits[child,]
      child_parent_d <- node_dist[parent,i]
      child_parent_w <- nexp(child_parent_d)
      
      if (is.null(pred) && is.null(total_w)){
        
        total_w <- total_w + child_parent_w
        pred <- child_trait * child_parent_w
        
      }else{
        
        pred <- pred + child_trait * child_parent_w
        total_w <- total_w + child_parent_w
        
      }
      
    }
    
  }
  
  if (is.null(pred)) return(NULL)
  
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

create_tree_and_traits <- function(fasta,ncores=1,verbose=FALSE){
  
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
  d <- dist.ml(alignment)
  if (verbose) cat('complete.\n')
  
  if (verbose) cat('Creating tree... ')
  tree <- nj(d) 
  if (verbose) cat('complete.\n')
  
  if (verbose) cat('Formatting tree and traits... ')
  tree <- multi2di(tree,random=FALSE)
  tree$node.label <- paste0('internal_node_',seq_along(unique(tree$edge[,1])))
  tree$edge.length[tree$edge.length < .0001] <- .0001
  
  img_ko <- read.delim(file.path(attr(fasta,'path'),'IMG_ko_counts.tab'),sep='\t')
  otu_img <- read.delim(file.path(attr(fasta,'path'),'gg_13_5_img.txt'),sep='\t')
  
  otu_img <- otu_img[!duplicated(otu_img$img_genome_id,fromLast=TRUE),]
  rownames(otu_img) <- otu_img$img_genome_id
  otu_img <- otu_img[otu_img$img_genome_id %in% img_ko$GenomeID,]
  otu_img <- otu_img[otu_img$X.gg_id %in% tree$tip.label,]
  
  tree_pruned <- drop.tip(tree,tip=tree$tip.label[!(tree$tip.label %in% otu_img$X.gg_id)])
  
  traits_pruned <- img_ko[img_ko$GenomeID %in% otu_img$img_genome_id,]
  rownames(traits_pruned) <- otu_img[as.character(traits_pruned$GenomeID),'X.gg_id']
  traits_pruned <- traits_pruned[,-1]
  
  out <- list(traits_pruned=traits_pruned,tree_pruned=tree_pruned,tree=tree)
  if (verbose) cat('complete.\n')
  
  return(out)
  
}

asr <- function(tree,traits){
  
  traits <- traits[tree$tip.label,]
  
  reconstructions <- apply(traits,2,ace,tree,type='continuous',method='pic') 
  
  asr_traits <- lapply(seq_along(reconstructions),function(x) reconstructions[[x]]$ace)
  names(asr_traits) <- names(reconstructions)
  
  asr_traits <- do.call(cbind,asr_traits)
  rownames(asr_traits) <- tree$node.label
  
  return(asr_traits)
  
}

predict_traits <- function(traits,traits_asr,tree){
  
  traits_asr <- readRDS('~/rcrust/asr.rds')
  traits <- traits[,colnames(traits_asr)]
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
tree_and_traits <- create_tree_and_traits(fasta,ncores=50,verbose=TRUE)
asr_traits <- asr(tree_and_traits$tree_pruned,tree_and_traits$traits_pruned)
reconstruction <- predict_traits(tree_and_traits$traits_pruned,asr_traits,tree_and_traits$tree)



# answer
final_picrust_table <- read.delim('~/rcrust/new_refs/ko_13_5_precalculated.tab',sep='\t')
rownames(final_picrust_table) <- final_picrust_table$X.OTU_IDs
final_picrust_table <- final_picrust_table[,-1]
cbind(rowSums(final_picrust_table[sort(rownames(final_picrust_table)[grepl('study',rownames(final_picrust_table))]),]),
      rowSums(reconstruction[sort(rownames(reconstruction)[grepl('rsv',rownames(reconstruction))]),]))


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

