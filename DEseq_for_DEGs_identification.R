rm(list=ls())
gc()
source("http://bioconductor.org/biocLite.R") 
biocLite("DESeq") 
library(DESeq)
count_table=read.table(file.choose(),sep="\t",header=TRUE,row.names=1)
count_table1<-apply(count_table, 1, function(hu) sum(hu<1)/length(hu)<=0.5)
count_table2<-count_table[count_table1,]
head(count_table2)
count_table2[,-1]<- lapply(lapply(count_table2[,-1],round),as.integer) #trandform the float into integer
m1 <- as.matrix(count_table2) #trandform the dataset into matrix
condition=factor(c(rep("Hyper",61),rep("Non_hyper",322))) # group the samples
storage.mode(m1) = "integer"
cds=newCountDataSet(m1,condition)
cds=estimateSizeFactors( cds )
cds=estimateDispersions( cds )
res=nbinomTest( cds, "Hyper", "Non_hyper" ) # Differential expression test
res_sig<-subset(res, padj<0.05&abs(log2FoldChange)>1)
write.table(res_sig,file="D:/CRC_DEseq_result.txt",quote=F,sep="\t",row.names=F)
