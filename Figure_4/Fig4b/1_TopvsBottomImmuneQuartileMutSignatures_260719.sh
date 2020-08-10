


################################################################
#     investigate the enrichment of mutational signatures in
#   samples with the highest vs lowest immune scores
#   perform wilcoxon's test on the signature weights
################################################################
cd ~/work/tumor_project/data/redo_samples/9_Strelka2_indels

#grab all signature scores
x=$(mktemp)
ls deconstructSigs_signatures_breast_only/SD*txt | while read l; do
  sample=$(basename $l | sed 's/.txt//')
  cat $l | awk '(NR%2==1){
    line1=line1"\t"$0
  }
  (NR%2==0){
    sample=$1
    $1=""
    line2=line2"\t"$0
  }
  END{
    print line1
    print line2
  }' | sed 's/ /\t/g' | \
  awk '(NR==1){for(i=1;i<=NF;i++) sig[i]=$i}
  (NR==2){for(i=1;i<=NF;i++) val[i]=$i}
  END{
    for(i=1;i<=NF;i++)
    print sig[i],val[i]
  }' > ${x}_${sample}
done

wc -l ${x}_* | awk '$1!=13'
#all have 13 lines, thus they all report the same number of data entries

#limit to those with >=15 SNVs
#and have both WES tumor and normal
awk 'NR>1 && $2=="." && $3=="."{print $1}' ../final_samples_090519.txt | sed 's/SD0855/SD0855_Rt/' | while read l; do
  if [ -s ${l}.filtered.vcf.gz ]; then
    echo -en $l"\t"$(zcat ${l}.filtered.vcf.gz | grep -v "#" | awk 'length($4)==1 && length($5)==1' | wc -l)"\n"
  fi
done | awk '$2>=15' > ${x}_totake
#506 samples

#create matrix of the samples vs signatures
(#print header
  l=SD0012
  cat ${x}_${l} | awk '{l=l"\t"$1} END{print l}'

  cut -f1 ${x}_totake | while read l; do
  echo $l $(cat ${x}_${l} | awk '{l=l"\t"$2} END{print l}')
done  ) > ${x}_grandtable

#take only samples for which ESTIMATE scores and signatures are available
#rename SD0855_Rt to SD0855
awk -vOFS="\t" '$1=="SD0855_Rt"{$1="SD0855"}{print}' ${x}_grandtable | \
  awk -vp=../../../annotations/immune_score/SJMC_immunescores.txt \
    'BEGIN{
      while((getline<p)>0)
        score[$1]=$8}
        (NR>1){print $1,score[$1]}' | awk 'NF==2 && $2!="NA"' > ${x}_samplesOK
#477 samples

#sort samples by ESTIMATE scores
#take only samples that are OK
#take top vs bottom quartile
quartile_size=$(cat ${x}_samplesOK | wc -l | awk '{print int($1/4+0.5)}')
#119
sort -k2,2gr ${x}_samplesOK | head -n $quartile_size > ${x}_top
sort -k2,2gr ${x}_samplesOK | tail -n $quartile_size > ${x}_bottom

#build the tables with the samples from the different categories
for y in top bottom; do
  cat ${x}_grandtable | \
    awk -vp=${x}_${y} 'BEGIN{
      while((getline<p)>0) k[$1]=1
    }
    ($1 in k) || NR==1' > ${x}_grandtable_${y}
done

wc -l ${x}_grandtable_*

echo $x
/tmp/tmp.lvKeuQOdMb





########
