


y="/tmp/tmp.CSssKCHzdS"

#calculate percentage
#Asians
sjmc=read.table(paste0(y,"_sjmc"))
smc=read.table(paste0(y,"_smc"))
tcga_asian=read.table(paste0(y,"_tcga_asian"))
nikzainal_korean=read.table(paste0(y,"_nikzainal_korean"))
tcga_caucasian=read.table(paste0(y,"_tcga_caucasian"))
nikzainal_nonkorean=read.table(paste0(y,"_nikzainal_nonkorean"))

sjmc_pct=sjmc$V2/sjmc$V3*100
smc_pct=smc$V2/smc$V3*100
nikzainal_korean_pct=nikzainal_korean$V2/nikzainal_korean$V3*100
tcga_asian_pct=tcga_asian$V2/tcga_asian$V3*100
nikzainal_nonkorean_pct=nikzainal_nonkorean$V2/nikzainal_nonkorean$V3*100
tcga_caucasian_pct=tcga_caucasian$V2/tcga_caucasian$V3*100

#create table
pct_all=rbind(sjmc_pct,smc_pct,tcga_asian_pct,nikzainal_korean_pct,
tcga_caucasian_pct,nikzainal_nonkorean_pct)
colnames(pct_all)=substr(sjmc$V1,11,20)

#plot all datasets

col_datasets=c(rgb(143,37,50,max=255),rgb(236,104,20,max=255),rgb(244,188,27,max=255),rgb(216,173,116,max=255),
rgb(116,117,162,max=255),rgb(175,186,216,max=255))

pdf(paste0("/home/mzabidi/all13sigs_alldatasets.pdf"),useDingbats=FALSE)

barstats=barplot(pct_all,col=col_datasets,xlab="signatures",ylab="% samples",
ylim=c(0,110),beside=TRUE,border=NA)

#draw legend
legend_text=c(paste0(rownames(pct_all)[1]," (",sjmc$V3[1],")"),
paste0(rownames(pct_all)[2]," (",smc$V3[1],")"),
paste0(rownames(pct_all)[3]," (",tcga_asian$V3[1],")"),
paste0(rownames(pct_all)[4]," (",nikzainal_korean$V3[1],")"),
paste0(rownames(pct_all)[5]," (",tcga_caucasian$V3[1],")"),
paste0(rownames(pct_all)[6]," (",nikzainal_nonkorean$V3[1],")"),
"P-value 2-sided","P-value g","P-value l")

legend("topright",legend_text,pch=15,bty="n",col=c(col_datasets,"white","white","white"))

Asian_pct=rbind(sjmc_pct,smc_pct,tcga_asian_pct,nikzainal_korean_pct)
Caucasian_pct=rbind(tcga_caucasian_pct,nikzainal_nonkorean_pct)

for(i in c(1:length(sjmc_pct)))
{
  p_res=round(t.test(Asian_pct[,i],Caucasian_pct[,i])$p.val,4)
  p_res_g=round(t.test(Asian_pct[,i],Caucasian_pct[,i],alternative="g")$p.val,4)
  p_res_l=round(t.test(Asian_pct[,i],Caucasian_pct[,i],alternative="l")$p.val,4)
  text(barstats[3,i],max(Asian_pct[,i],Caucasian_pct[,i])+5,
  paste0(p_res,"\n",p_res_g,"\n",p_res_l),cex=0.8)
}

dev.off()





###########
