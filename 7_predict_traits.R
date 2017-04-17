nexp <- function(d) exp(1)^(-d)

method <- 'asr_and_weighting'
weighting <- 'exponential'

asr_picrust <- read.delim('~/rcrust/asr/KEGG_asr_counts.tab',sep='\t')

reconstructed_trait_table <- readRDS('~/rcrust/asr.rds')
tree_full <- read.tree('~/rcrust/out/KEGG/reference_tree.newick')
observed_trait_table <- read.delim('out/KEGG/trait_table.tab',sep='\t')
rownames(observed_trait_table) <- observed_trait_table$GenomeID
observed_trait_table <- observed_trait_table[,-1]

traits <- rbind(observed_trait_table,reconstructed_trait_table)





test <- read.delim('~/rcrust/new_refs/ko_13_5_precalculated.tab',sep='\t')
test_steps <- read.delim('~/MiscOut/testTable.tab',sep='\t')



preorder <- reorder.phylo(tree_full,order='cladewise')
