


R
library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )
holder="/tmp/tmp.icm6EN21f8"

#n breaks
my.breaks <- c(0,0.5,1,2,3)
#n-1 colors
my.colors<-c("gray","#E11E26","blue","lightyellow3")
my.legend=c("amplification","deletion","NOT DONE")

pdf(paste0("/home/mzabidi/IntClust_4_1.65cutoff.pdf"))

e<-read.table(paste0(holder,"_FINAL"))

title="mutations in each sample and types"
heatmap.2(as.matrix(e),
  main=title,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace="none",
  density.info="none",cexCol=0.2,cexRow=1,
  ,margins = c(5,5),key=TRUE,symkey=TRUE,
  breaks=my.breaks,scale="none",col=my.colors)

legend("bottomleft",my.legend,col=my.colors[-1],pch=15,bty="n")

dev.off()













########
