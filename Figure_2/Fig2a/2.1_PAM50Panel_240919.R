

library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )

#n breaks
my.breaks <- c(0,0.9,1,2,3,4,5)
#n-1 colors
my.colors<-c("gray","royalblue","plum2","red3","olivedrab3","lightslategray")
my.legend=c("NA","Luminal A","Luminal B","Her2-enriched","Basal-like","Normal-like")

pdf(paste0("/home/mzabidi/PAM50.pdf"))

e<-read.table("/tmp/tmp.mHFcgLP3Md_ready_BIG",header=T)

nmatrix=cbind(e$PAM50,e$PAM50_DUMMY)
colnames(nmatrix)=c("PAM50","PAM50dummy")
row.names(nmatrix)=row.names(e)

title="PAM50"
heatmap.2(t(nmatrix),
  main=title,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace="none",
  density.info="none",cexCol=0.2,cexRow=1,
  ,margins = c(5,5),key=TRUE,symkey=TRUE,
  breaks=my.breaks,scale="none",col=my.colors)

legend("bottomleft",my.legend,col=my.colors,pch=15,bty="n")


dev.off()







########
