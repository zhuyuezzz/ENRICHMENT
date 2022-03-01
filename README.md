# ENRICHMENT
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")


setwd("C:\\Users\\lenovo\\Desktop\\GSE25097 - TRPM8\\15.symbo2id")    

library("org.Hs.eg.db")                                           
rt=read.table("diff.txt",sep="\t",check.names=F,header=T)         
rt=rt[,c(1,2)]
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)       
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)       



#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05  

qvalueFilter=0.05                   

setwd("D:\\biowolf\\geoSurvival\\16.GO")          
zz=read.table("id.txt",sep="\t",header=T,check.names=F)           
zz=zz[is.na(zz[,"entrezID"])==F,]                                 
gene=zz$entrezID
geneFC=zz$logFC
names(geneFC)=gene

colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}


kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)                
showNum=30
if(nrow(GO)<30){
	showNum=nrow(GO)
}

pdf(file="barplot.pdf",width = 11,height = 7)
barplot(kk, drop = TRUE, showCategory =showNum,color = colorSel,font.size = 10)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=150))
dev.off()		#关闭pdf
pdf(file="bubble.pdf",width = 11,height = 7)
dotplot(kk,showCategory = showNum, orderBy = "GeneRatio", color = colorSel,font.size = 10)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))+ theme_bw() +
  theme(panel.grid=element_blank())



#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05                
qvalueFilter=0.05                


zy=read.table("id.txt",sep="\t",header=T,check.names=F)            
zy=zy[is.na(zy[,"entrezID"])==F,]                                  
gene=zy$entrezID
geneFC=zy$logFC
names(geneFC)=gene

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


ss <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)   
KEGG=as.data.frame(ss)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(zy$id[match(strsplit(x,"/")[[1]],as.character(zy$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)                          
read.table("KEGG.txt",sep="\t",header = T,check.names=F)
showNum=30
if(nrow(KEGG)<30){
  showNum=nrow(KEGG)
}


pdf(file="barplot.pdf",width = 8,height = 7)
barplot(ss, drop = TRUE, showCategory = showNum, color = colorSel,font.size = 10)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60)) + theme_bw() +
  theme(panel.grid=element_blank())
dev.off()


pdf(file="bubble.pdf",width = 8,height = 7)
dotplot(ss, showCategory = showNum, orderBy = "GeneRatio",color = colorSel,font.size=8)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

