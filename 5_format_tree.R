library(ape)

tree_pruned_picrust <- read.tree('out/KEGG/pruned_tree.newick')
trait_table <- read.delim('out/KEGG/trait_table.tab',sep='\t')

tree_full <- read.tree('out/KEGG/reference_tree.newick')
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

saveRDS(list(traits=traits_pruned,tree=tree_pruned),'~/rcrust/pruned.rds')


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