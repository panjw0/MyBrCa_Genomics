

library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )

library(gplots)
library(MASS)
library(RColorBrewer)

my.breaks <- c(seq(-0.3,0.3,0.005))
#my.breaks <- c(seq(-0.5,0.5,0.005))
#n-1 colors
my.colors<-colorpanel(length(my.breaks)-1,rgb(49,105,154,25,maxColorValue=255),
    rgb(255,255,255,25,maxColorValue=255),rgb(164,25,32,25,maxColorValue=255))          #x-1 colors

#paste the results of echo ${x}_grandtable_immune here:
d=read.table("/tmp/tmp.OiREvF4APY_grandtable_immune")

#remove NA data points
take_data_IMPRES=which(!is.na(d$IMPRES))
take_data_IFNg_GSVA=which(!is.na(d$IFNg_GSVA))
take_data_Bindea_GSVA=which(!is.na(d$Bindea_GSVA))
take_data_ESTIMATE=which(!is.na(d$ESTIMATE))

way="spearman"
i=1
w=cor(d[take_data_IMPRES,i],d$IMPRES[take_data_IMPRES],method=way)
x=cor(d[take_data_IFNg_GSVA,i],d$IFNg_GSVA[take_data_IFNg_GSVA],method=way)
y=cor(d[take_data_Bindea_GSVA,i],d$Bindea_GSVA[take_data_Bindea_GSVA],method=way)
z=cor(d[take_data_ESTIMATE,i],d$ESTIMATE[take_data_ESTIMATE],method=way)
for(i in seq(2,13))
{
  w=rbind(w,cor(d[take_data_IMPRES,i],d$IMPRES[take_data_IMPRES],method=way))
  x=rbind(x,cor(d[take_data_IFNg_GSVA,i],d$IFNg_GSVA[take_data_IFNg_GSVA],method=way))
  y=rbind(y,cor(d[take_data_Bindea_GSVA,i],d$Bindea_GSVA[take_data_Bindea_GSVA],method=way))
  z=rbind(z,cor(d[take_data_ESTIMATE,i],d$ESTIMATE[take_data_ESTIMATE],method=way))
}

f=cbind(w,x,y,z)
colnames(f)<-c("IMPRES","IFNg_GSVA","Bindea_GSVA","ESTIMATE")
rownames(f)<-c(colnames(d)[1:13])

pdf(paste0("/home/mzabidi/immune_score_vs_signature_weights_",way,"_0.3max.pdf"))

title=paste0(way," corr sig weights\nvs immune score\n")
heatmap.2(as.matrix(f),
main=title,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace="none",
density.info="none",margins = c(9,9),cexCol=1,cexRow=1,
key=TRUE,symkey=TRUE,
breaks=my.breaks,scale="none",col=my.colors,na.color="gray")

dev.off()


######################################################
#draw scatter plot to see how the data look like
######################################################
#do for ESTIMATE
my.colors=c(rainbow(16))

y=d$ESTIMATE[take_data_ESTIMATE]
for(i in seq(1,16))
{
  x=d[take_data_ESTIMATE,i]

  pdf(paste0("/home/mzabidi/ESTIMATE_vs_sig_",colnames(d)[i],".pdf"),useDingbats=F)
  plot(x,y,xlab=paste0(colnames(d)[i]," weight"),ylab="ESTIMATE score",
  main=paste0(colnames(d)[i]," vs ESTIMATE score"),col=my.colors[i],frame.plot=FALSE,pch=19)

  pcc=round(cor(x,y,method="pearson"),4)
  pcc_pval=round(cor.test(x,y,method="pearson")$p.val,6)
  scc=round(cor(x,y,method="spearman"),4)
  scc_pval=round(cor.test(x,y,method="spearman")$p.val,6)
  kcc=round(cor(x,y,method="kendall"),4)
  kcc_pval=round(cor.test(x,y,method="kendall")$p.val,6)

  legend("topright",paste0("PCC=",pcc,"\t",pcc_pval,
  "\nSCC=",scc,"\t",scc_pval,
  "\nKCC=",kcc,"\t",kcc_pval),bty="n")

  dev.off()

}




########
