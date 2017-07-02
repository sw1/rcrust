library(ape)
library(phangorn)

options(digits=10) 

nexp <- function(x) exp(1)^(-x) # -1x

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
    
    child <- tree$tip.label[i]
    
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
  
  pred <- unlist(pred/total_w)
  
  return(pred)
  
}


# prepare traits
trait_asr_table <- read.delim('/data/sw1/rcrust/asr/KEGG_asr_counts.tab',sep='\t',stringsAsFactors=FALSE)
rownames(trait_asr_table) <- trait_asr_table[,1]
trait_asr_table <- trait_asr_table[,-1]

trait_table <- read.delim('/data/sw1/rcrust/out/KEGG/trait_table.tab',sep='\t',stringsAsFactors=FALSE)
rownames(trait_table) <- trait_table[,1]
trait_table <- trait_table[,-1]
trait_table <- trait_table[,colnames(trait_asr_table)]

traits <- rbind(trait_asr_table,trait_table)



# get nodes to predict
pitree <- read.tree('~/rcrust/tree_before_predictnode.tree')

# picrust has internal_node_0 for root, probably remove and keep as root for mine
pitree$node.label[which(pitree$node.label == "'internal_node_0'")] <- 'root'


node_dist <- dist.nodes(pitree)

profvis({

recons <- matrix(0.0,length(pitree$tip.label),ncol(traits),dimnames=list(pitree$tip.label,colnames(traits)))
for (i in 1:5){ #seq_along(pitree$tip.label)){
  
  node_to_predict <- pitree$tip.label[i] # pitree$tip.label == nodes_to_predict
  
  cat(sprintf('%s: Tip %s... ',i,node_to_predict))
  
  # reconstruction
  ancestors <- Ancestors(pitree,i,type='all')
  mrra <- NULL
  for (j in ancestors){
    ancestor <- pitree$node.label[j - length(pitree$tip.label)]
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
    recons[node_to_predict,] <- unlist(traits[node_to_predict,])
  }else{
    recons[node_to_predict,] <- pred
  }
  
}


})





# 
# Ancestors(pitree,10,type='all')
# Children(pitree,Ancestors(pitree,10,type='all'))
# 
# lapply(Ancestors(pitree,10,type='all'),function(x) pitree$edge[pitree$edge[,1] %in% x,])
# lapply(Ancestors(pitree,10,type='all'),function(x) pitree$edge.length[pitree$edge[,1] %in% x])
