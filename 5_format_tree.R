library(ape)
library(readr)

tree_pruned <- read.tree('out/KEGG/pruned_tree.newick')
trait_table <- read.delim('out/KEGG/trait_table.tab',sep='\t')

tree_full <- read.tree('out/KEGG/reference_tree.newick')
keggGenome_counts <- read.delim('picrust_starting_files/IMG_ko_counts.tab',sep='\t')
otu_to_keggGenome <- read.delim('gg_13_5_img.txt',sep='\t')



# picrust style

trait_table_fields <- keggGenome_counts

trait_to_tree_mapping <- vector(mode='list')
for (i in seq_len(nrow((otu_to_keggGenome)))){
  key <- as.character(otu_to_keggGenome$img_genome_id[i])
  value <- as.character(otu_to_keggGenome$`X.gg_id`[i])
  trait_to_tree_mapping[key] <- value
}

remapped <- vector(mode='list')
for (i in seq_len(nrow(trait_table_fields))){
  field <- as.character(unlist(trait_table_fields[i,1]))
  remapped <- c(remapped,trait_to_tree_mapping[[field]])
}
remapped <- unlist(remapped)

parsed <- remapped[remapped %in% tree_full$tip.label]
length(parsed)



# R style

mapping_unique <- otu_to_keggGenome[!duplicated(otu_to_keggGenome$img_genome_id,fromLast=TRUE),]
otu_intersect <- mapping_unique[mapping_unique$img_genome_id %in% keggGenome_counts$GenomeID,'X.gg_id']
parsed <- otu_intersect[otu_intersect %in% tree_full$tip.label]
