library(ape)

asr_picrust <- read.delim('~/rcrust/asr/KEGG_asr_counts.tab',sep='\t')

pruned <- readRDS('~/rcrust/pruned.rds')
tree <- pruned$tree
data <- pruned$traits

data_ordered <- data[tree$tip.label,-1]

reconstructions <- apply(data_ordered,2,ace,tree,type="continuous",method='pic') #default to pic

#pull out only the ace node predictions
just_ace <- lapply(seq_along(reconstructions),function(x) reconstructions[[x]]$ace)
names(just_ace) <- names(reconstructions)

#reformat the list into a matrix
just_ace_matrix <- do.call(cbind,just_ace)

#relabel the node names (ones created internally by ape) with the actual node labels in the tree
# just_ace_matrix <- cbind(tree$node.label,just_ace_matrix)
# colnames(just_ace_matrix)[1]<-'nodes'
# out_matrix <- as.data.frame(just_ace_matrix,check.names=FALSE)
# write.table(out_matrix,file=count_out_file,row.names=FALSE,quote=FALSE, sep="\t")
rownames(just_ace_matrix) <- tree$node.label # picrust names the 0 node "root"

saveRDS(just_ace_matrix,'~/rcrust/asr.rds')




# 
# 
# #extract 95% CI info
# ci <- lapply(seq_along(reconstructions),
#              function(x) paste(round(reconstructions[[x]]$CI95[,1],digits=4),
#                                round(reconstructions[[x]]$CI95[,2],digits=4),sep="|"))
# names(ci)<-names(reconstructions)                                                                                
# ci_matrix<-do.call(cbind,ci)
# ci_matrix<-cbind(tree$node.label,ci_matrix)
# 
# if(asr_method=="ML" || asr_method=="REML"){
#   #extract information about brownian motion parameter
#   sigma<-lapply(1:length(reconstructions),function(x) paste(round(reconstructions[[x]]$sigma2[1],digits=4),round(reconstructions[[x]]$sigma2[2],digits=4),sep="|"))
#   names(sigma)<-names(reconstructions)                                                                                
#   sigma_matrix<-do.call(cbind,sigma)
#   sigma_matrix<-cbind(c('sigma'),sigma_matrix)
#   
#   #add it to the ci matrix
#   ci_matrix<-rbind(ci_matrix,sigma_matrix)
#   
#   #get loglik value
#   if(asr_method=="ML"){
#     loglik<-lapply(1:length(reconstructions),function(x) round(reconstructions[[x]]$loglik,digits=4))                                                                                                                           
#   }else{
#     loglik<-lapply(1:length(reconstructions),function(x) round(reconstructions[[x]]$resloglik,digits=4))                                                                                                                                                                                        
#   }
#   
#   names(loglik)<-names(reconstructions)
#   loglik_matrix<-do.call(cbind,loglik)
#   loglik_matrix<-cbind(c('loglik'),loglik_matrix)
#   ci_matrix<-rbind(ci_matrix,loglik_matrix)
#   
# }
# 
# #just set the column name for the row names to someting arbtrary like 'nodes'
# colnames(ci_matrix)[1]<-'nodes'
# 
# #output the data to file
# out_matrix<-data.frame(ci_matrix,check.names=FALSE)                                                                                                                 
# write.table(out_matrix,file=ci_out_file,row.names=FALSE,quote=FALSE, sep="\t")