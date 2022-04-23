library(DESeq2)
setwd("C:\\Users\\Harleen\\Downloads\\Task1\\Task1")
countdata<-read.csv("Count_data.csv",header = TRUE, )
countdata
metadata<-read.csv("Metadata.csv",header = TRUE, )
metadata
counts<-countdata[,2:8]
row.names(counts) <- countdata[,1]
dd <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design= ~ condition)
dim(dd)
dd <- dds[ rowSums(counts(dds)) > 10, ]
dd


ddse <- estimateSizeFactors(dd)
se <- SummarizedExperiment(log2(counts(ddse, normalized=TRUE) + 1),
                           colData=colData(ddse))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
plotPCA( DESeqTransform( se ) )
dds<-DESeq(dd)
dds
re<-results(dds,alpha = 0.05)
re
summary(re)
deseq_results<-results(dds, contrast=c("condition", "mock", "dht"))
write.csv(deseq_results,"deseq_results.csv")
