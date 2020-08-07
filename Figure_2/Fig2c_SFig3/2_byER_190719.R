

##################################################################
# barplot of mutational prevalence from the different datasets
#               according to ER status
##################################################################
#Asians
SJMC_file="/home/mzabidi/work/tumor_project/data/redo_samples/9_Strelka2_indels/ERstatus_mutation_MyBRCA.txt"
ZhengyanKan2018_file="/home/mzabidi/work/tumor_project/data/published_data/ZhengyanKan2018/ERstatus_mutation_Korean.txt"
TCGA_Asian_file="/home/mzabidi/work/tumor_project/data/published_data/brca_tcga/ERstatus_mutation_TCGA_Asians.txt"
NikZainalkorean_file="/home/mzabidi/work/tumor_project/data/published_data/NikZainal_291018/nextera_regions/ERstatus_mutation_korean_NikZainal2016.txt"
#Caucasians
TCGA_Caucasian_file="/home/mzabidi/work/tumor_project/data/published_data/brca_tcga/ERstatus_mutation_TCGA_Caucasians.txt"
METABRIC_file="/home/mzabidi/work/tumor_project/data/published_data/brca_metabric/ERstatus_mutation_METABRIC.txt"
NikZainalcaucasian_file="/home/mzabidi/work/tumor_project/data/published_data/NikZainal_291018/nextera_regions/ERstatus_mutation_nonkorean_NikZainal2016.txt"

#read file
#Asians
SJMC_data=read.table(SJMC_file)
ZhengyanKan2018_data=read.table(ZhengyanKan2018_file)
TCGA_Asian_data=read.table(TCGA_Asian_file)
NikZainalkorean_data=read.table(NikZainalkorean_file)
#non-Asians
TCGA_Caucasian_data=read.table(TCGA_Caucasian_file)
METABRIC_data=read.table(METABRIC_file)
NikZainalcaucasian_data=read.table(NikZainalcaucasian_file)

#calculate percentage
#Asians
SJMC_pct=SJMC_data$all_subtypes[-1]/SJMC_data$all_subtypes[1]*100
ZhengyanKan2018_pct=ZhengyanKan2018_data$all_subtypes[-1]/ZhengyanKan2018_data$all_subtypes[1]*100
TCGA_Asian_pct=TCGA_Asian_data$all_subtypes[-1]/TCGA_Asian_data$all_subtypes[1]*100
NikZainalkorean_pct=NikZainalkorean_data$all_subtypes[-1]/NikZainalkorean_data$all_subtypes[1]*100

#non-Asians
TCGA_Caucasian_pct=TCGA_Caucasian_data$all_subtypes[-1]/TCGA_Caucasian_data$all_subtypes[1]*100
METABRIC_pct=METABRIC_data$all_subtype[-1]/METABRIC_data$all_subtype[1]*100
NikZainalcaucasian_pct=NikZainalcaucasian_data$all_subtypes[-1]/NikZainalcaucasian_data$all_subtypes[1]*100

#create table
pct_all=rbind(SJMC_pct,ZhengyanKan2018_pct,TCGA_Asian_pct,NikZainalkorean_pct,
METABRIC_pct,TCGA_Caucasian_pct,NikZainalcaucasian_pct)

colnames(pct_all)=c(rownames(SJMC_data)[-1])

col_datasets=c(rgb(143,37,50,max=255),rgb(236,104,20,max=255),rgb(244,188,27,max=255),rgb(216,173,116,max=255),
rgb(139,140,233,max=255),rgb(116,117,162,max=255),rgb(175,186,216,max=255))

pdf(paste0("/home/mzabidi/topgenes_ERstatus_alldatasets.pdf"),useDingbats=FALSE)
par(mfrow=c(2,1))

for(h in c(2:3))
{
  subtype=colnames(SJMC_data)[h]
  SJMC_pct=SJMC_data[-1,h]/SJMC_data[1,h]*100
  ZhengyanKan2018_pct=ZhengyanKan2018_data[-1,h]/ZhengyanKan2018_data[1,h]*100
  TCGA_Asian_pct=TCGA_Asian_data[-1,h]/TCGA_Asian_data[1,h]*100
  NikZainalkorean_pct=NikZainalkorean_data[-1,h]/NikZainalkorean_data[1,h]*100

  METABRIC_pct=METABRIC_data[-1,h]/METABRIC_data[1,h]*100
  TCGA_Caucasian_pct=TCGA_Caucasian_data[-1,h]/TCGA_Caucasian_data[1,h]*100
  NikZainalcaucasian_pct=NikZainalcaucasian_data[-1,h]/NikZainalcaucasian_data[1,h]*100

  pct_subtype=rbind(SJMC_pct,ZhengyanKan2018_pct,TCGA_Asian_pct,NikZainalkorean_pct,
  METABRIC_pct,TCGA_Caucasian_pct,NikZainalcaucasian_pct)
  colnames(pct_subtype)=c(rownames(SJMC_data)[-1])

  barstats=barplot(pct_subtype,col=col_datasets,xlab="Genes",
    ylab=paste0("% ",subtype),border=NA,beside=TRUE,
    ylim=c(0,110),cex.names=0.8,
    main=paste0("Asian vs Caucasian comparison\n ER ",subtype,"-ive"))

  #plot legend
  legend_text=c(paste0(rownames(pct_all)[1]," (",SJMC_data[1,h],")"),
    paste0(rownames(pct_all)[2]," (",ZhengyanKan2018_data[1,h],")"),
    paste0(rownames(pct_all)[3]," (",TCGA_Asian_data[1,h],")"),
    paste0(rownames(pct_all)[4]," (",NikZainalkorean_data[1,h],")"),
    paste0(rownames(pct_all)[5]," (",METABRIC_data[1,h],")"),
    paste0(rownames(pct_all)[6]," (",TCGA_Caucasian_data[1,h],")"),
    paste0(rownames(pct_all)[7]," (",NikZainalcaucasian_data[1,h],")"),
    "P-val 2-sided t-test")

  legend("topright",legend_text,pch=15,bty="n",col=c(col_datasets,"white"))

  #write t-test p-values
  Asian_pct=rbind(SJMC_pct,ZhengyanKan2018_pct,TCGA_Asian_pct,NikZainalkorean_pct)
  Caucasian_pct=rbind(METABRIC_pct,TCGA_Caucasian_pct,NikZainalcaucasian_pct)
  colnames(Asian_pct)=c(rownames(SJMC_data)[-1])
  colnames(Caucasian_pct)=c(rownames(SJMC_data)[-1])

  for(i in c(1:length(rownames(SJMC_data)[-1])))
  {
    p_res=round(t.test(Asian_pct[,i],Caucasian_pct[,i])$p.val,6)
    text(barstats[3,i],max(Asian_pct[,i],Caucasian_pct[,i])+5,
    paste0("P=",p_res),cex=0.8)
  }

  #"

}

dev.off()



#also adjust ylim and make it bigger
h=2     #ERnegative
for(ylim_height in c(100,10,5))
h=3     #ERpositive
{
  #pdf(paste0("/home/mzabidi/R_output/topgenes_ERneg_alldatasets_ylim",ylim_height,".pdf"),useDingbats=FALSE)
  pdf(paste0("/home/mzabidi/R_output/topgenes_ERpos_alldatasets_ylim",ylim_height,".pdf"),useDingbats=FALSE)

  subtype=colnames(SJMC_data)[h]
  SJMC_pct=SJMC_data[-1,h]/SJMC_data[1,h]*100
  ZhengyanKan2018_pct=ZhengyanKan2018_data[-1,h]/ZhengyanKan2018_data[1,h]*100
  TCGA_Asian_pct=TCGA_Asian_data[-1,h]/TCGA_Asian_data[1,h]*100
  NikZainalkorean_pct=NikZainalkorean_data[-1,h]/NikZainalkorean_data[1,h]*100

  METABRIC_pct=METABRIC_data[-1,h]/METABRIC_data[1,h]*100
  TCGA_Caucasian_pct=TCGA_Caucasian_data[-1,h]/TCGA_Caucasian_data[1,h]*100
  NikZainalcaucasian_pct=NikZainalcaucasian_data[-1,h]/NikZainalcaucasian_data[1,h]*100

  pct_subtype=rbind(SJMC_pct,ZhengyanKan2018_pct,TCGA_Asian_pct,NikZainalkorean_pct,
  METABRIC_pct,TCGA_Caucasian_pct,NikZainalcaucasian_pct)
  colnames(pct_subtype)=c(rownames(SJMC_data)[-1])

  barstats=barplot(pct_subtype,col=col_datasets,xlab="Genes",
    ylab=paste0("% ",subtype),border=NA,beside=TRUE,
    ylim=c(0,ylim_height),cex.names=0.8,
    main=paste0("Asian vs Caucasian comparison\n ER ",subtype,"-ive"))

  #plot legend
  legend_text=c(paste0(rownames(pct_all)[1]," (",SJMC_data[1,h],")"),
    paste0(rownames(pct_all)[2]," (",ZhengyanKan2018_data[1,h],")"),
    paste0(rownames(pct_all)[3]," (",TCGA_Asian_data[1,h],")"),
    paste0(rownames(pct_all)[4]," (",NikZainalkorean_data[1,h],")"),
    paste0(rownames(pct_all)[5]," (",METABRIC_data[1,h],")"),
    paste0(rownames(pct_all)[6]," (",TCGA_Caucasian_data[1,h],")"),
    paste0(rownames(pct_all)[7]," (",NikZainalcaucasian_data[1,h],")"),
    "P-val 2-sided t-test")

  legend("topright",legend_text,pch=15,bty="n",col=c(col_datasets,"white"))

  #write t-test p-values
  Asian_pct=rbind(SJMC_pct,ZhengyanKan2018_pct,TCGA_Asian_pct,NikZainalkorean_pct)
  Caucasian_pct=rbind(METABRIC_pct,TCGA_Caucasian_pct,NikZainalcaucasian_pct)
  colnames(Asian_pct)=c(rownames(SJMC_data)[-1])
  colnames(Caucasian_pct)=c(rownames(SJMC_data)[-1])

  for(i in c(1:length(rownames(SJMC_data)[-1])))
  {
    p_res=round(t.test(Asian_pct[,i],Caucasian_pct[,i])$p.val,6)
    text(barstats[3,i],max(Asian_pct[,i],Caucasian_pct[,i])+5,
    paste0("P=",p_res),cex=0.8)
  }

  dev.off()

}



###########
