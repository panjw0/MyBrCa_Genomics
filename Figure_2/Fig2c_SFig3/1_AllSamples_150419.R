


##################################################################
#barplot of mutational prevalence from the different datasets
#           all samples regardless of subtypes
##################################################################

#Asians
SJMC_file="/people/mzabidi/work/tumor_project/data/redo_samples/9_Strelka2_indels/subtype_mutation_MyBRCA.txt"
ZhengyanKan2018_file="/people/mzabidi/work/tumor_project/data/published_data/ZhengyanKan2018/subtype_mutation_Korean.txt"
NikZainal2016_Korean_file="/people/mzabidi/work/tumor_project/data/published_data/NikZainal_291018/nextera_regions/subtype_mutation_korean_NikZainal2016.txt"
TCGA_Asian_file="/people/mzabidi/work/tumor_project/data/published_data/brca_tcga/subtype_mutation_TCGA_Asians.txt"
#Caucasians
TCGA_Caucasian_file="/people/mzabidi/work/tumor_project/data/published_data/brca_tcga/subtype_mutation_TCGA_Caucasians.txt"
METABRIC_file="/people/mzabidi/work/tumor_project/data/published_data/brca_metabric/subtype_mutation_METABRIC.txt"
NikZainal2016_NonKorean_file="/people/mzabidi/work/tumor_project/data/published_data/NikZainal_291018/nextera_regions/subtype_mutation_nonkorean_NikZainal2016.txt"


#read file
#Asians
SJMC_data=read.table(SJMC_file)
ZhengyanKan2018_data=read.table(ZhengyanKan2018_file)
NikZainal2016_Korean_data=read.table(NikZainal2016_Korean_file)
TCGA_Asian_data=read.table(TCGA_Asian_file)
#non-Asians
TCGA_Caucasian_data=read.table(TCGA_Caucasian_file)
METABRIC_data=read.table(METABRIC_file)
NikZainal2016_NonKorean_data=read.table(NikZainal2016_NonKorean_file)


#calculate percentage
#Asians
SJMC_pct=SJMC_data$all_subtypes[-1]/SJMC_data$all_subtypes[1]*100
ZhengyanKan2018_pct=ZhengyanKan2018_data$all_subtypes[-1]/ZhengyanKan2018_data$all_subtypes[1]*100
NikZainal2016_Korean_pct=NikZainal2016_Korean_data$all_subtypes[-1]/NikZainal2016_Korean_data$all_subtypes[1]*100
TCGA_Asian_pct=TCGA_Asian_data$all_subtypes[-1]/TCGA_Asian_data$all_subtypes[1]*100
#non-Asians
TCGA_Caucasian_pct=TCGA_Caucasian_data$all_subtypes[-1]/TCGA_Caucasian_data$all_subtypes[1]*100
METABRIC_pct=METABRIC_data$all_subtype[-1]/METABRIC_data$all_subtype[1]*100
NikZainal2016_NonKorean_pct=NikZainal2016_NonKorean_data$all_subtypes[-1]/NikZainal2016_NonKorean_data$all_subtypes[1]*100

#create table
pct_all=rbind(SJMC_pct,ZhengyanKan2018_pct,TCGA_Asian_pct,NikZainal2016_Korean_pct,
METABRIC_pct,TCGA_Caucasian_pct,NikZainal2016_NonKorean_pct)
colnames(pct_all)=c(rownames(SJMC_data)[-1])

#plot all datasets
col_datasets=c(rgb(143,37,50,max=255),rgb(236,104,20,max=255),rgb(244,188,27,max=255),rgb(216,173,116,max=255),
rgb(139,140,233,max=255),rgb(116,117,162,max=255),rgb(175,186,216,max=255))

pdf(paste0("/home/mzabidi/top10genes_alldatasets.pdf"),useDingbats=FALSE)
barstats=barplot(pct_all,col=col_datasets,xlab="Genes",ylab="% all samples",
beside=TRUE,border=NA,
ylim=c(0,80),cex.names=0.6,
main="Asian vs Caucasian comparison\nall samples")

#draw legend
legend_text=c(paste0(rownames(pct_all)[1]," (",SJMC_data$all_subtypes[1],")"),
paste0(rownames(pct_all)[2]," (",ZhengyanKan2018_data$all_subtypes[1],")"),
paste0(rownames(pct_all)[3]," (",TCGA_Asian_data$all_subtypes[1],")"),
paste0(rownames(pct_all)[4]," (",NikZainal2016_Korean_data$all_subtypes[1],")"),
paste0(rownames(pct_all)[5]," (",METABRIC_data$all_subtype[1],")"),
paste0(rownames(pct_all)[6]," (",TCGA_Caucasian_data$all_subtypes[1],")"),
paste0(rownames(pct_all)[7]," (",NikZainal2016_NonKorean_data$all_subtypes[1],")"),
"P-val 2-sided t-test","1-sided greater","1-sided less")

legend("topright",legend_text,pch=15,bty="n",col=c(col_datasets,"white","white","white"))

Asian_pct=rbind(SJMC_pct,ZhengyanKan2018_pct,TCGA_Asian_pct,NikZainal2016_Korean_pct)
Caucasian_pct=rbind(METABRIC_pct,TCGA_Caucasian_pct,NikZainal2016_NonKorean_pct)
colnames(Asian_pct)=c(rownames(SJMC_data)[-1])
colnames(Caucasian_pct)=c(rownames(SJMC_data)[-1])

for(i in c(1:length(rownames(SJMC_data)[-1])))
{
  p_res=round(t.test(Asian_pct[,i],Caucasian_pct[,i])$p.val,4)
  p_res_g=round(t.test(Asian_pct[,i],Caucasian_pct[,i],alternative="g")$p.val,4)
  p_res_l=round(t.test(Asian_pct[,i],Caucasian_pct[,i],alternative="l")$p.val,4)
  text(barstats[3,i],max(Asian_pct[,i],Caucasian_pct[,i])+5,
  paste0(p_res,"\n",p_res_g,"\n",p_res_l),cex=0.8)
}

dev.off()


###########
