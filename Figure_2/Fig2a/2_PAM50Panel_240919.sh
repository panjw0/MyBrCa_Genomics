

#################################################################
#        plot PAM50 classification bar panels
#  use the sorting of samples as per IntClust classification
#################################################################
cd ~/work/tumor_project/data/redo_samples/9_Strelka2_indels/annotated_MAF

#grab sample names, separated by IntClusts
samples=$(mktemp)
awk '$0~/^#/{for(i=2;i<=NF;i++) print $i}' compiled_clean.maf > ${samples}_all
cat ../../../../annotations/Study_metadata_8.5.19.csv | \
  awk -vFS="," -vp=$samples 'NR>1{print $1 > p"_"$18}'
#546 samples

#the top 10 genes
#also add other interesting genes (11) at the back
genes=$(mktemp)
echo -en "TP53\nPIK3CA\nGATA3\nMAP3K1\nKMT2C\nPTEN\nCBFB\nCDH1\nAKT1\nNF1\n" > $genes
echo -en "TBX3\nARID1A\nRB1\nMYH8\nUSP9X\nNOTCH1\nNOTCH2\nRUNX1\nERBB3\nNCOR1\nPALB2\n" >> $genes

#create matrix of samples with their mutation types
#separate according to IntClusts
matrix=$(mktemp)
for IntClust in $(awk -vFS="," 'NR>1{print $18}' ../../../../annotations/Study_metadata_8.5.19.csv | sort | uniq); do
  echo $IntClust
  echo -en "\t"$(cat $genes | sed 's/\n/\t/')"\n" > ${matrix}_header_${IntClust}

  for curr_sample in $(cat ${samples}_${IntClust}); do
    awk '$1!~"#"' compiled_clean.maf | \
      awk -vp=$curr_sample -vFS="\t" -vq=$genes \
      'BEGIN{
        #grab genes to consider
        while((getline<q)>0)
          genes[++totgenes]=$1
        close(q)
        #empty the gene statistics
        for(i=1;i<=totgenes;i++)
          gene_stat[i]=0
      }
      ($16==p){
        #here, change the annotation to Multi_Hit if >1 and different mutations for the gene in the same sample
        for(i=1;i<=totgenes;i++)
          if($1==genes[i])
            if(!gene_stat[i] || (gene_stat[i] && gene_stat[i]==$9)) {gene_stat[i]=$9}
          else {gene_stat[i]="Multi_Hit"}
      }
      END{
        line=p
        for(i=1;i<=totgenes;i++)
          line=line"\t"gene_stat[i]
        print line
      }'
  done > ${matrix}_raw_${IntClust}
done

#prepare term dictionary, matching the type of mutation to internal code
dict=$(mktemp)
echo -en "Frame_Shift_Ins\t1
Frame_Shift_Del\t1
Splice_Site\t1
Nonsense_Mutation\t1
Start_Codon_Del\t1
Start_Codon_SNP\t1
In_Frame_Ins\t2
In_Frame_Del\t2
Missense_Mutation\t2
Multi_Hit\t3
0\t0\n" > $dict

#translate the terms into the code above
for IntClust in $(awk -vFS="," 'NR>1{print $18}' ../../../../annotations/Study_metadata_8.5.19.csv | sort | uniq); do
  cat ${matrix}_raw_${IntClust} | \
    awk -vp=$dict 'BEGIN{
      while((getline<p)>0)
      code[$1]=$2
    }
    {
      line=$1
      for(i=2;i<=NF;i++)
      line=line"\t"code[$i]
      print line
    }' > ${matrix}_trans_${IntClust}
done

#create sorting template for the samples
for IntClust in $(awk -vFS="," 'NR>1{print $18}' ../../../../annotations/Study_metadata_8.5.19.csv | sort | uniq); do
  cat ${matrix}_raw_${IntClust} | \
    awk -vOFS="\t" '{
      line=$1
      for(i=2;i<=NF;i++)
      line=line"\t"($i?1:0)
      print line
    }' | sort -k2,2r -k3,3r -k4,4r -k5,5r -k6,6r -k7,7r -k8,8r -k9,9r -k10,10r -k11,11r > ${matrix}_sort_template_${IntClust}
done

#make the final matrix to input in R
for IntClust in $(awk -vFS="," 'NR>1{print $18}' ../../../../annotations/Study_metadata_8.5.19.csv | sort | uniq); do
  ( cat ${matrix}_header_${IntClust}
    cat ${matrix}_sort_template_${IntClust} | \
    awk -vp=${matrix}_trans_${IntClust} 'BEGIN{
      while((getline<p)>0)
      line[$1]=$0
    }
    {print line[$1]}' ) > ${matrix}_ready_${IntClust}
done

#grab PAM50 and create the final matrix to input in R
for IntClust in $(awk -vFS="," 'NR>1{print $18}' ../../../../annotations/Study_metadata_8.5.19.csv | sort | uniq); do
  ( echo -en "\tPAM50\tPAM50_DUMMY\n"
    cat ${matrix}_sort_template_${IntClust} | \
      awk -vp=../../../../annotations/immune_score/SJMC_immunescores.txt 'BEGIN{
        while((getline<p)>0)
          PAM50[$1]=$3
      }
      {
        c="NA"
        if(PAM50[$1]=="LumA") c=1
        if(PAM50[$1]=="LumB") c=2
        if(PAM50[$1]=="Her2") c=3
        if(PAM50[$1]=="Basal") c=4
        if(PAM50[$1]=="Normal") c=5

        print $1,c,"NA"}' ) > ${matrix}_ready_${IntClust}
done

#concatenate the matrices one after another
#put separators
IntClust=1
separator_cols=$(cat ${matrix}_ready_${IntClust} | awk 'NR==1{for(i=1;i<=NF;i++) $i="NA"; print}')
( cat ${matrix}_ready_${IntClust}
  echo -en "ENDOF_"$IntClust"_1\t"$separator_cols"\n"
  echo -en "ENDOF_"$IntClust"_2\t"$separator_cols"\n"  ) > ${matrix}_ready_BIG
for IntClust in 2 3 4- 4+ 5 6 7 8 9 10 NA; do
  awk 'NR>1' ${matrix}_ready_${IntClust}
  echo -en "ENDOF_"$IntClust"_1\t"$separator_cols"\n"
  echo -en "ENDOF_"$IntClust"_2\t"$separator_cols"\n"
done >> ${matrix}_ready_BIG


echo $matrix



R
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
