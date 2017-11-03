library(ShortRead)

ref <- readFasta("~/rcrust/gg_13_5.fasta.gz")
img_ko <- read.table(file = "~/rcrust/picrust_starting_files/IMG_ko_counts.tab", sep = "\t", row.names = 1, header = TRUE, comment.char = '')
id_img <- read.table(file = "~/rcrust/gg_13_5_img.txt", sep = "\t", header = TRUE, comment.char = '')

img_ko_sub <- img_ko[rownames(img_ko) %in% id_img$img_genome_id, ]
id_img_sub <- id_img[id_img$img_genome_id %in% rownames(img_ko_sub), ]

ids <- id(ref)
ref_sub <- ref[as.character(ids) %in% id_img_sub$X.gg_id, ]
ids_db <- as.character(id(ref_sub))
seqs_db <- as.character(sread(ref_sub))

seqtab.nochim <- t(readRDS('~/rcrust/seqtab.rds'))
seqs_study <- rownames(seqtab.nochim)
ids_study <- paste("study", 1:nrow(seqtab.nochim), sep = "_")
# merge db and study seqs
ids_out <- c(ids_db, ids_study)
seqs_out <- c(seqs_db, seqs_study)
fasta <- ShortRead(sread = DNAStringSet(seqs_out), id = BStringSet(ids_out))
writeFasta(fasta, file = "~/rcrust/gg_13_5_study_db.fasta")
# (optional) output gg id-to-IMG mapping info; modify this to filter later for predicted traits (predict_traits.py -l )
write.table(as.character(id(fasta)), file = "~/rcrust/traits_sample_filter.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# 2) output seq variant count data as biom;
# note: use Joey's biom latest dev version; library(devtools); install_github("joey711/biom")

saveRDS(fasta,'~/rcrust/fasta.rds')
