logFoldChange=1
adjustP=0.05

library(limma)
setwd("GSE12021")
rt=read.table("normalize.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)


modType=c(rep("RA",12),rep("Ctrl",9),rep("RA",13),rep("Ctrl",10))
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("Ctrl","RA")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(RA-Ctrl,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=T,row.names=T)

diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=T,row.names=T)
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=T,row.names=T)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=T,row.names=T)


hmExp=rt[as.vector(diffSig[,1]),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)
