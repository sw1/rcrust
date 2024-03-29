library(ape)


# picrust to compare
tree_pruned_picrust <- read.tree('out/KEGG/pruned_tree.newick')
trait_table <- read.delim('out/KEGG/trait_table.tab',sep='\t')






# R
# tree_full <- read.tree('~/rcrust/out/KEGG/reference_tree.newick')
tree_full <- read.tree('~/rcrust/tree.tree')
#tree_full <- root(tree_full,tree_full$tip.label[1],resolve.root=TRUE)
tree_full <- multi2di(tree_full,random=FALSE)
tree_full$node.label <- paste0('internal_node_',seq_along(unique(tree_full$edge[,1])))
tree_full$edge.length[tree_full$edge.length < .0001] <- .0001
# tree_full <- root(reorder.phylo(tree_full,order='postorder'),length(tree_full$tip.label),resolve.root=TRUE)

keggGenome_counts <- read.delim('picrust_starting_files/IMG_ko_counts.tab',sep='\t')
otu_to_keggGenome <- read.delim('gg_13_5_img.txt',sep='\t')

# R style

mapping_unique <- otu_to_keggGenome[!duplicated(otu_to_keggGenome$img_genome_id,fromLast=TRUE),]
rownames(mapping_unique) <- mapping_unique$img_genome_id
otu_intersect <- mapping_unique[mapping_unique$img_genome_id %in% keggGenome_counts$GenomeID,]
parsed <- otu_intersect[otu_intersect$X.gg_id %in% tree_full$tip.label,]

tree_pruned <- drop.tip(tree_full,tip=tree_full$tip.label[!(tree_full$tip.label %in% parsed$X.gg_id)])

traits_pruned <- keggGenome_counts[keggGenome_counts$GenomeID %in% parsed$img_genome_id,]
rownames(traits_pruned) <- parsed[as.character(traits_pruned$GenomeID),'X.gg_id']

saveRDS(list(traits=traits_pruned,tree=tree_pruned,tree_full=tree_full),'~/rcrust/pruned.rds')
write.tree(tree_pruned,'~/rcrust/pruned.tree')
write.tree(tree_full,'~/rcrust/tree_full.tree')

trait_table_r <- traits_pruned
trait_table_r[,1] <- rownames(traits_pruned)
rownames(trait_table_r) <- NULL
write.table(trait_table_r,'~/rcrust/trait_table_r.tab',sep='\t')



# # picrust style
# 
# trait_table_fields <- keggGenome_counts
# 
# trait_to_tree_mapping <- vector(mode='list')
# for (i in seq_len(nrow((otu_to_keggGenome)))){
#   key <- as.character(otu_to_keggGenome$img_genome_id[i])
#   value <- as.character(otu_to_keggGenome$`X.gg_id`[i])
#   trait_to_tree_mapping[key] <- value
# }
# 
# remapped <- vector(mode='list')
# for (i in seq_len(nrow(trait_table_fields))){
#   field <- as.character(unlist(trait_table_fields[i,1]))
#   remapped <- c(remapped,trait_to_tree_mapping[[field]])
# }
# remapped <- unlist(remapped)
# 
# parsed <- remapped[remapped %in% tree_full$tip.label]
# length(parsed)