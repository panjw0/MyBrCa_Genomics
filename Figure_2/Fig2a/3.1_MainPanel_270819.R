


R
library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )
holder="/tmp/tmp.KdJ7rIf4EP"

#n breaks
my.breaks <- c(0,0.5,1,2,3)
#n-1 colors
my.colors<-c("gray","#E11E26","#35A047","black")
my.legend=c("FSindels_nonsense_splice_startloss","IFindels_missense","Multi-hit")

#pdf(paste0("/home/mzabidi/IntClust.pdf"))
pdf(paste0("/home/mzabidi/IntClust_Pereirra_arrangement.pdf"))

e<-read.table(paste0(holder,"_ready_BIG"))

title="mutations in each sample and types"
heatmap.2(as.matrix(t(e)),
  main=title,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace="none",
  density.info="none",cexCol=0.2,cexRow=1,
  ,margins = c(5,5),key=TRUE,symkey=TRUE,
  breaks=my.breaks,scale="none",col=my.colors)

legend("bottomleft",my.legend,col=my.colors[-1],pch=15,bty="n")

dev.off()




########
