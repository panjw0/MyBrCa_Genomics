

library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )

holder="/tmp/tmp.lvKeuQOdMb"
top_sigs=read.table(paste0(holder,"_grandtable_top"))
bottom_sigs=read.table(paste0(holder,"_grandtable_bottom"))

pdf(paste0("/home/mzabidi/p-vals_wilcoxon_top_vs_bottom.pdf"))

par(mfrow=c(1,2))

#higher in the lower immune score
i=1
p_val=wilcox.test(top_sigs[,i],bottom_sigs[,i],alternative="less")$p.val
for(i in c(seq(2,13)))
{p_val=c(p_val,wilcox.test(top_sigs[,i],bottom_sigs[,i],alternative="less")$p.val)}

#adjust p-value using bonferroni
p_adjusted=p.adjust(p_val, method = "bonferroni", n = length(p_val))

barplot(rev(-log(p_adjusted,10)),col=rev(rainbow(13)),border=NA,names.arg=substr(rev(colnames(top_sigs)),11,15),
ylab="signatures",xlab="-log10(P-value Wilcoxon's test)",
main="-log10 P-value(Bonferroni-corr) 1-sided Wilcoxon's(higher for lower ESTIMATE)
of signatures in samples with top 100 vs bottom 100 ESTIMATE scores",
space=0.05,horiz=T,xlim=c(0,6))

abline(v=-log(cutoff,10),col="grey",lwd=3,lty=2)
text(-log(cutoff,10),"P=0.05")


#higher in the higher immune score
cutoff=0.05
i=1
p_val=wilcox.test(top_sigs[,i],bottom_sigs[,i],alternative="greater")$p.val
for(i in c(seq(2,13)))
{p_val=c(p_val,wilcox.test(top_sigs[,i],bottom_sigs[,i],alternative="greater")$p.val)}

#adjust p-value using bonferroni
p_adjusted=p.adjust(p_val, method = "bonferroni", n = length(p_val))

barplot(rev(-log(p_adjusted,10)),col=rev(rainbow(13)),border=NA,names.arg=substr(rev(colnames(top_sigs)),11,15),
ylab="signatures",xlab="-log10(P-value Wilcoxon's test)",
main="-log10 P-value(Bonferroni-corr) 1-sided Wilcoxon's(higher for higher ESTIMATE)
of signatures in samples with top 100 vs bottom 100 ESTIMATE scores",
space=0.05,horiz=T,xlim=c(0,6))

abline(v=-log(cutoff,10),col="grey",lwd=3,lty=2)
text(-log(cutoff,10),"P=0.05")

dev.off()


#later manually arrange the signatures according to Fig.4A





########
