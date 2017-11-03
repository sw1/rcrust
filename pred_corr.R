# Prediction corr

abundance <- readRDS('~/Documents/rcrust/rsv.rds')
colnames(abundance) <- paste('rsv',1:ncol(abundance),sep='_')

abundance_norm <- round(t(t(abundance)/reconstruction_16s[colnames(abundance),]))

fxns <- reconstruction_ko[colnames(abundance_norm),]
ids <- intersect(colnames(abundance_norm),rownames(fxns))

fxns <- fxns[ids,]
abundance_norm <- abundance_norm[,ids]

predictions <- round(abundance_norm %*% fxns)
rownames(predictions) <- gsub('\\_[0-9]$','',rownames(predictions))

truth <- read.delim('~/Downloads/kos.tsv.gz',sep='\t')
rownames(truth) <- truth$X..Gene.Family
truth <- truth[,-1]
truth <- truth[!grepl('\\|',rownames(truth)),]
truth <- truth[-c(1:2),]

mapping <- read.csv('~/Documents/rcrust/wgs_amp_mapping.txt',stringsAsFactors=FALSE)


pred_ss <- predictions[rownames(predictions) %in% mapping$amp_err,]
rownames(pred_ss) <- mapping$wgs_id[match(rownames(pred_ss),mapping$amp_err)]
truth_ss <- t(truth)
rownames(truth_ss) <- mapping$wgs_id[match(gsub('^(SRR[0-9]+)_Abundance.RPKs$','\\1',rownames(truth_ss)),mapping$wgs_srr)]
truth_ss <- truth_ss[rownames(pred_ss),]

truth_ss <- truth_ss/rowSums(truth_ss)
pred_ss <- pred_ss/rowSums(pred_ss)

x <- truth_ss[2,intersect(colnames(truth_ss),colnames(pred_ss))]
y <- pred_ss[2,intersect(colnames(truth_ss),colnames(pred_ss))]
cor(x,y)
qplot(x,y,geom='point',alpha=.4) + geom_smooth(method='lm',color='red') + scale_x_log10() + scale_y_log10()


all_kos <- unique(c(colnames(pred_ss),colnames(t_ss)))

ko_missing_pred <- all_kos[!(all_kos %in% colnames(pred_ss))]
ko_missing_pred_mat <- matrix(0.0,nrow(pred_ss),length(ko_missing_pred),dimnames=list(rownames(pred_ss),ko_missing_pred))
full_pred_ss <- cbind(pred_ss,ko_missing_pred_mat)

ko_missing_t <- all_kos[!(all_kos %in% colnames(t_ss))]
ko_missing_t_mat <- matrix(0.0,nrow(t_ss),length(ko_missing_t),dimnames=list(rownames(t_ss),ko_missing_t))
full_t_ss <- cbind(t_ss,ko_missing_t_mat)
full_t_ss <- full_t_ss[rownames(full_pred_ss),colnames(full_pred_ss)]

# full correlation
cor(c(full_t_ss),c(full_pred_ss))

# Absence, presents corr
cor(as.integer(c(full_t_ss) > 0),as.integer(c(full_pred_ss) > 0),method='spearman')

# overlapping entries correlation
overlap <- intersect(names(which(colSums(full_t_ss) > 0)),names(which(colSums(full_pred_ss) > 0)))
cor(c(full_t_ss[,overlap]),c(full_pred_ss[,overlap]))


plot(log10(c(full_t_ss)+.00001),log10(c(full_pred_ss)+.00001))
