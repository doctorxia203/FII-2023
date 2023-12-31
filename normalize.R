library(sva)
library(limma)
setwd("GSE12021")
rt=read.table("merge.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
batchType=c(rep(1,48),rep(2,6))
modType=c(rep("Ctrl",24),rep("RIF",24),rep("Ctrl",3),rep("RIF",3))
mod = model.matrix(~as.factor(modType))
outTab=ComBat(data, batchType, mod, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="normalize.txt",sep="\t",quote=F,col.names=F)