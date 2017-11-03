library(ape)
library(phangorn)
library(DECIPHER)

align_seqs <- '/usr/local/bin/align_seqs.py'

fasta <- readRDS('~/Documents/rcrust/fasta.rds')

### do alignment
seqs <- sread(fasta)
names(seqs) <- as.character(id(fasta))

alignment <- AlignSeqs(sread(fasta),processors=NULL)
names(alignment) <- as.character(id(fasta))

alignment_phydat <- phyDat(as(alignment,'matrix'),type='DNA')

### compute pairwise distances 
dm <- dist.ml(alignment_phydat)

### perform neighbor-joining tree estimation
tree <- NJ(dm) 

### save tree
ape::write.tree(tree,'~/rcrust/tree.tree')



# QIIME style
system2(align_seqs,
        args=c('-e','90',
               '-p','0.1',
               '-i','~/rcrust/gg_13_5_study_db.fasta',
               '-o','~/rcrust/gg_13_5_study_db_aligned'))

system('/usr/local/bin/FastTree -nt -gamma -fastest -no2nd -spr 4 ~/rcrust/gg_13_5_study_db_aligned/gg_13_5_study_db_aligned.fasta > ~/rcrust/qiime_tree.tree')
        