library(ShortRead)
library(ape)
library(phangorn)
library(DECIPHER)


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



abundance <- readRDS('~/Documents/rcrust/seqtab.rds')
path <- '~/Documents/rcrust/reference_files'

# change dist.ml
# change nj

fasta <- create_fasta(abundance,path,FALSE,verbose=TRUE)
saveRDS(fasta,'/home/dag332/rcrust/alignment_large.rds')