



R
library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )
holder="/tmp/tmp.wTKz7fah0p"

#n breaks
my.breaks <- c(0,0.5,2.5,2.8)
#n-1 colors
my.colors<-c("gray","black","lightyellow3")
my.legend=c("NEGATIVE","POSITIVE","NOTKNOWN")

pdf(paste0("/home/mzabidi/HRpanel.pdf"))

e<-read.table(paste0(holder,"_ready_BIG"))

title="PR, ER, HER2, and HER2 final status"
heatmap.2(as.matrix(t(e)),
  main=title,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace="none",
  density.info="none",cexCol=0.2,cexRow=1,
  ,margins = c(5,5),key=TRUE,symkey=TRUE,
  breaks=my.breaks,scale="none",col=my.colors)

legend("bottomleft",my.legend,col=my.colors,pch=15,bty="n")

dev.off()







########
