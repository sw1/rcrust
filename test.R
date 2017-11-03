truth <- read.delim('~/Documents/rcrust/reference_files/Picrust/ko_tail_final.tab', sep='\t')
recon <- readRDS('~/Documents/rcrust/testing_workflow/reconstruction_ko.rds')

ko_ids <- intersect(colnames(recon),colnames(truth)[-c(1,ncol(truth))])

recon <- recon[,ko_ids]
recon <- cbind('#OTU_IDs'=rownames(recon),recon,'metadata_NSTI'='0.0')
recon <- rbind(recon,as.matrix(truth[,colnames(recon)]))

write_delim(data.frame(recon,stringsAsFactors=FALSE),
            '~/Documents/rcrust/reference_files/Picrust/ko_recon.tab.gz',col_names=TRUE,delim='\t')


######################

abundance <- readRDS('~/Documents/rcrust/rsv.rds')
colnames(abundance) <- paste('rsv',1:ncol(abundance),sep='_')

abundance_norm <- round(t(t(abundance)/reconstruction_16s[colnames(abundance),]))

fxns <- reconstruction_ko[colnames(abundance_norm),]
ids <- intersect(colnames(abundance_norm),rownames(fxns))

fxns <- fxns[ids,]
abundance_norm <- abundance_norm[,ids]

predictions <- round(abundance_norm %*% fxns)
predictions_norm <- predictions/rowSums(predictions)
