dir_name <- '~/Documents/rcrust/testing_workflow/dist_jc69_tree_nj'
trees <- list.files(dir_name,full.names=TRUE)


mat <- matrix(0,choose(length(trees),2),6,dimnames=list(NULL,c('i','j','sym','br','path','quad')))
df <- as.data.frame(matrix(NA,length(trees),length(trees),dimnames=list(1:length(trees),1:length(trees))))
row <- 1
for (i in seq(1,length(trees)-1)){
  for (j in seq(i+1,length(trees))){
    
    tree1 <- read.tree(trees[i])
    tree2 <- read.tree(trees[j])
    
    d <- treedist(tree1,tree2)
    
    mat[row,] <- c(i,j,d)
    df[[i,j]] <- data.frame(d)
    
    row <- row + 1
    
  }
}

#expand.grid(d=c('ml','raw','N','TS','TV', 'JC69', 'K80', 'F81','K81','F84','BH87','T92','TN93','GG95','logdet','paralin','indel','indelblock'),
#            tree=c('nj','bionj','fastme.ols','fastme.bal'))
