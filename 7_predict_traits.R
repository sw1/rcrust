library(ape)
library(phangorn)

options(digits=10) 

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


# answer
final_picrust_table <- read.delim('~/rcrust/new_refs/ko_13_5_precalculated.tab',sep='\t')
rownames(final_picrust_table) <- final_picrust_table$X.OTU_IDs
final_picrust_table <- final_picrust_table[,-1]





# prepare traits (from r 6_asr.R)
trait_asr_table <- readRDS('~/rcrust/asr.rds')

trait_table <- read.delim('/data/sw1/rcrust/out/KEGG/trait_table.tab',sep='\t',stringsAsFactors=FALSE)
rownames(trait_table) <- trait_table[,1]
trait_table <- trait_table[,-1]
trait_table <- trait_table[,colnames(trait_asr_table)]

prepared_data <- readRDS('~/rcrust/pruned.rds')
trait_table <- prepared_data$traits[,-1]
trait_table <- trait_table[,colnames(trait_asr_table)]

traits <- as.matrix(rbind(trait_asr_table,trait_table))

pitree <- prepared_data$tree_full




# # fix answer
# final_picrust_table <- as.matrix(final_picrust_table)
# rownames(final_picrust_table)[grepl('study',rownames(final_picrust_table))] <- paste0("'",rownames(final_picrust_table)[grepl('study',rownames(final_picrust_table))],"'")
# 
# # prepare traits (from picrust format_tree_and_trait_table)
# trait_asr_table <- read.delim('/data/sw1/rcrust/asr/KEGG_asr_counts.tab',sep='\t',stringsAsFactors=FALSE)
# rownames(trait_asr_table) <- trait_asr_table[,1]
# trait_asr_table <- trait_asr_table[,-1]
# 
# trait_table <- read.delim('/data/sw1/rcrust/out/KEGG/trait_table.tab',sep='\t',stringsAsFactors=FALSE)
# rownames(trait_table) <- trait_table[,1]
# trait_table <- trait_table[,-1]
# trait_table <- trait_table[,colnames(trait_asr_table)]
# 
# traits <- as.matrix(rbind(trait_asr_table,trait_table))
# 
# # this is the preprocessed tree from picrust right before predict nodes
# # get nodes to predict
# pitree <- read.tree('~/rcrust/tree_before_predictnode.tree')
# # # picrust has internal_node_0 for root, probably remove and keep as root for mine
# pitree$node.label[which(pitree$node.label == "'internal_node_0'")] <- 'root'






node_dist <- dist.nodes(pitree)

recons <- matrix(0.0,length(pitree$tip.label),ncol(traits),dimnames=list(pitree$tip.label,colnames(traits)))
for (i in seq_along(pitree$tip.label)){
  
  # node_to_predict <- pitree$tip.label[i] # pitree$tip.label == nodes_to_predict
  node_to_predict <- get_node_label(pitree,i)
  
  cat(sprintf('%s: Tip %s... ',i,node_to_predict))
  
  # reconstruction
  ancestors <- Ancestors(pitree,i,type='all')
  mrra <- NULL
  for (j in ancestors){
    # ancestor <- pitree$node.label[j - length(pitree$tip.label)]
    ancestor <- get_node_label(pitree,j)
    if (ancestor %in% rownames(traits)){
      mrra <- ancestor
      mrra_idx <- j
      break
    }
  }
  
  cat(sprintf('predicting with mrra %s\n',mrra))
  pred <- round(predict_tip(pitree,i,mrra,mrra_idx,traits,node_dist))
  
  
  # update reconstruction table
  if (node_to_predict %in% rownames(traits)){
    recons[node_to_predict,] <- traits[node_to_predict,]
  }else{
    recons[node_to_predict,] <- pred
  }
  
}




cbind(rowSums(recons[rownames(recons)[grepl('study',rownames(recons))],]),
      rowSums(final_picrust_table[rownames(recons)[grepl('study',rownames(recons))],]))




# Ancestors(pitree,10,type='all')
# Children(pitree,Ancestors(pitree,10,type='all'))
# 
# lapply(Ancestors(pitree,10,type='all'),function(x) pitree$edge[pitree$edge[,1] %in% x,])
# lapply(Ancestors(pitree,10,type='all'),function(x) pitree$edge.length[pitree$edge[,1] %in% x])
