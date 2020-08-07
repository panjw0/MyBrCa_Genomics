

#############################################################
#           plot hormone receptor panels
#############################################################
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
  ( echo -en "\tPR\tER\tHER2\tHER2final\n"
    cat ${matrix}_sort_template_${IntClust} | \
      awk -vp=../../../../annotations/Study_metadata_8.5.19.csv 'BEGIN{
        FS=","
        while((getline<p)>0)
        {
          PR[$1]=$7
          ER[$1]=$8
          HER2[$1]=$9
          HER2final[$1]=$11
        }
        FS="\t"
      }
      {print $1,PR[$1],ER[$1],HER2[$1],HER2final[$1]}' | \
        awk -vOFS="\t" '{
          for(i=2;i<=5;i++)
          {
            #0 is negative, 3 is NOT KNOWN, 1 is positive
            if($i=="N") $i="0"
            if($i=="NA") $i="3"
            if($i=="P" || $i=="1+" || $i=="2+" || $i=="3+") $i="1"
            if($i!="0" && $i!="NA" && $i!="1") $i="3"      #CHECK FOR REMAINING STUPIDITY
          }
          print
        }' ) > ${matrix}_ready_${IntClust}
done



#concatenate the matrices one after another
#put a separator
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
/tmp/tmp.wTKz7fah0p





########
