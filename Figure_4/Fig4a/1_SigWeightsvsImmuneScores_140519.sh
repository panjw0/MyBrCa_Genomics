

##########################################################################
#   heatmap of mutational signature weights vs immune score for MyBrCa
##########################################################################
cd ~/work/tumor_project/data/redo_samples/9_Strelka2_indels

#grab signature scores
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

#check data integrity
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
    if [ -s ${x}_${l} ]; then
      echo $l $(cat ${x}_${l} | awk '{l=l"\t"$2} END{print l}')
    fi
done  ) > ${x}_grandtable

#grab immune scores
awk -vOFS="\t" '$1=="SD0855_Rt"{$1="SD0855"}{print}' ${x}_grandtable | \
  awk -vOFS="\t" -vp=../../../annotations/immune_score/SJMC_immunescores.txt \
  'BEGIN{
    while((getline<p)>0)
    {
      IMPRES[$1]=$5
      IFNg_GSVA[$1]=$6
      Bindea_GSVA[$1]=$7
      ESTIMATE[$1]=$8
    }
  }
  (NR==1){print $0,"IMPRES","IFNg_GSVA","Bindea_GSVA","ESTIMATE"}
  (NR>1){print $0,IMPRES[$1]?IMPRES[$1]:"NA",IFNg_GSVA[$1]?IFNg_GSVA[$1]:"NA",Bindea_GSVA[$1]?Bindea_GSVA[$1]:"NA",ESTIMATE[$1]?ESTIMATE[$1]:"NA"}' | \
  grep -v NA > ${x}_grandtable_immune
#477 samples survived

echo ${x}_grandtable_immune
/tmp/tmp.OiREvF4APY



########
