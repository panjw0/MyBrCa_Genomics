


##########################################################################
#   heatmap of mutational signature weights vs immune score for TCGA
##########################################################################
cd ~/work/tumor_project/data/published_data/brca_tcga/vcfs

library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )

library(gplots)
library(MASS)
library(RColorBrewer)

#my.breaks <- c(seq(-0.3,0.3,0.005))
max=0.5
my.breaks <- c(seq(-0.5,0.5,0.005))


#n-1 colors
my.colors<-colorpanel(length(my.breaks)-1,rgb(49,105,154,25,maxColorValue=255),
    rgb(255,255,255,25,maxColorValue=255),rgb(164,25,32,25,maxColorValue=255))          #x-1 colors


for(race in c("Caucasian","Asian"))
{
  d=read.table(paste0("/tmp/tmp.Hrk6XFm5Ba_grandtable_immune_",race))

  #remove NA data points
  take_data_IMPRES=which(!is.na(d$IMPRES))
  take_data_eIFNgamma=which(!is.na(d$eIFNgamma))
  take_data_GSVABindea=which(!is.na(d$GSVABindea))
  take_data_ESTIMATE=which(!is.na(d$ESTIMATE))

  for(way in c("pearson","spearman","kendall"))
  {
    i=1
    w=cor(d[take_data_IMPRES,i],d$IMPRES[take_data_IMPRES],method=way)
    x=cor(d[take_data_eIFNgamma,i],d$eIFNgamma[take_data_eIFNgamma],method=way)
    y=cor(d[take_data_GSVABindea,i],d$GSVABindea[take_data_GSVABindea],method=way)
    z=cor(d[take_data_ESTIMATE,i],d$ESTIMATE[take_data_ESTIMATE],method=way)

    for(i in seq(2,13))
    {
      w=rbind(w,cor(d[take_data_IMPRES,i],d$IMPRES[take_data_IMPRES],method=way))
      x=rbind(x,cor(d[take_data_eIFNgamma,i],d$eIFNgamma[take_data_eIFNgamma],method=way))
      y=rbind(y,cor(d[take_data_GSVABindea,i],d$GSVABindea[take_data_GSVABindea],method=way))
      z=rbind(z,cor(d[take_data_ESTIMATE,i],d$ESTIMATE[take_data_ESTIMATE],method=way))
    }

    f=cbind(w,x,y,z)
    colnames(f)<-c("IMPRES","eIFNgamma","GSVABindea","ESTIMATE")
    rownames(f)<-c(colnames(d)[1:13])

    pdf(paste0("/home/mzabidi/",race,"_immune_score_vs_signature_weights_",way,"_",max,".pdf"))

    title=paste0(way," corr sig weights vs immune score\nTCGA ",race,"\n")
    heatmap.2(as.matrix(f),
    main=title,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace="none",
    density.info="none",margins = c(9,9),cexCol=1,cexRow=1,
    key=TRUE,symkey=TRUE,
    breaks=my.breaks,scale="none",col=my.colors,na.color="gray")

    dev.off()

  }

}


#####################


#draw scatter plot to see how the data look like
my.colors=c(rainbow(13))

for(race in c("Caucasian","Asian"))
{
  d=read.table(paste0("/tmp/tmp.z5DFKRsdXU_grandtable_immune_",race))

  #remove NA data points
  take_data_ESTIMATE=which(!is.na(d$ESTIMATE))

  y=d$ESTIMATE[take_data_ESTIMATE]
  for(i in seq(1,13))
  {
    x=d[take_data_ESTIMATE,i]

    pdf(paste0("/home/mzabidi/",race,"SGVA_vs_sig_",colnames(d)[i],".pdf"),useDingbats=F)
    plot(x,y,xlab=paste0(colnames(d)[i]," weight"),ylab="ESTIMATE score",
    main=paste0(colnames(d)[i]," vs ESTIMATE score\nTCGA ",race),col=my.colors[i],frame.plot=FALSE,pch=19)

    pcc=round(cor(x,y,method="pearson"),4)
    scc=round(cor(x,y,method="spearman"),4)
    kcc=round(cor(x,y,method="kendall"),4)

    legend("topright",paste0("PCC=",pcc,"\nSCC=",scc,"\nKCC=",kcc),bty="n")

    dev.off()
  }
}






########
