


######################################################################
#         plot a heatmap of amplification and deletion
#  use sample arrangement as per oncoplot for short mutations
######################################################################

cd ~/work/tumor_project/data/redo_samples/9_Strelka2_indels/annotated_MAF

############################################
#     grab sample sorting scheme
# here, use the sorting scheme of Fig2a
# hence have to re-run the code for Fig2a
# to rebuild the sorting scheme correctly
############################################
#grab sample names, separated by IntClusts
samples=$(mktemp)
awk '$0~/^#/{for(i=2;i<=NF;i++) print $i}' compiled_clean.maf > ${samples}_all
cat ../../../../annotations/Study_metadata_8.5.19.csv | \
  awk -vFS="," -vp=$samples 'NR>1{print $1 > p"_"$18}'
#546 samples

#the top 10 genes
genes=$(mktemp)
echo -en "TP53\nPIK3CA\nGATA3\nMAP3K1\nKMT2C\nPTEN\nCBFB\nCDH1\nAKT1\nNF1\n" > $genes
echo -en "TBX3\nARID1A\nRB1\nMYH8\nUSP9X\nNOTCH1\nNOTCH2\nRUNX1\nERBB3\nNCOR1\nPALB2" >> $genes

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

#make the final matrix (this what was input into R for the short mutations heatmap)
for IntClust in $(awk -vFS="," 'NR>1{print $18}' ../../../../annotations/Study_metadata_8.5.19.csv | sort | uniq); do
  ( cat ${matrix}_header_${IntClust}
    cat ${matrix}_sort_template_${IntClust} | \
    awk -vp=${matrix}_trans_${IntClust} 'BEGIN{
      while((getline<p)>0)
      line[$1]=$0
    }
    {print line[$1]}' ) > ${matrix}_ready_${IntClust}
done



##########################################################
#      plot amplification and deletion of genes
##########################################################

cna=$(mktemp)
genes=$(mktemp)
bed=$(mktemp)

#gene list
for gene in ERBB2 CCND1 MYC RPS6KB1 ZNF703 PTRH2 APPBP2 RSF1 INTS4 PPM1D PIK3CA TBL1XR1 \
  EMSY GATA3 KRAS CDKN1B AKT1 RB1 CDH1 TP53 CBFB MAP2K4 PPP2R2A; do
    echo $gene
  done > $genes

#determine amplification/deletion status
for IntClust in 1 2 3 4- 4+ {5..10} NA; do
  for sample in $(awk 'NR>1{print $1}' ${matrix}_ready_${IntClust}); do

    if [ -s ../../../copy_number/final_samples_segments/to_use/${sample}.txt ]; then

      echo $sample
      cat ../../../copy_number/final_samples_segments/to_use/${sample}.txt | \
        awk -vOFS="\t" 'NR>1{print $2,$3,$4,$5,".","."}' > ${bed}_${sample}.bed

      for gene in $(cat $genes); do

        line=$(zcat ../../../../annotations/UCSC_hg19_uniq.bed.gz | awk -vp=$gene '$4==p')
        chr=$(echo $line | awk '{print $1}' | sed 's/chr//')
        start=$(echo $line | awk '{print $2}')
        end=$(echo $line | awk '{print $3}')

        echo -en $chr"\t"$start"\t"$end"\n" > ${bed}_${gene}.bed

        #intersect, and take average value of CNA
        value=$(bedtools intersect -a  ${bed}_${sample}.bed -b ${bed}_${gene}.bed | awk '{t+=$4} END{print t/NR}')
          echo -en $gene"\t"$value"\n"
      done > ${cna}_${sample}

    else

      #if no sWGS done
      echo $sample "NO DATA"
      for gene in $(cat $genes); do
        echo -en $gene"\tNOSWGS\n"
      done > ${cna}_${sample}

    fi

  done
done


#transform into amplification or deletion or no-sWGS flags
for IntClust in {1..3} 4- 4+ {5..10} NA; do
  for sample in $(awk 'NR>1{print $1}' ${matrix}_ready_${IntClust}); do

    cat ${cna}_${sample} | \
      awk -vOFS="\t" -vp=$sample -vuppercut=3 -vundercut=1.65 \
        'BEGIN{
        print p
      }{
      CNA=0
      if(2^$2>uppercut) CNA=1
      if(2^$2<undercut) CNA=2
      if($2=="NOSWGS") CNA=3
      print CNA}' > ${cna}_${sample}_transformed

  done
done

for IntClust in {1..3} 4- 4+ {5..10} NA; do
  #create matrix
  pastelist=$(awk 'NR>1{print $1}' ${matrix}_ready_${IntClust} | while read sample; do
    echo ${cna}_${sample}_transformed
  done)
  paste $pastelist > ${matrix}_ready_${IntClust}_CNA

  #create filler columns (2 columns)
  cat ${matrix}_ready_${IntClust}_CNA | \
    awk -vOFS="\t" -vp=$IntClust \
      'NR==1{print "DUMMY_IntClust"p"_1","DUMMY_IntClust"p"_2"}
        NR>1{print "NA","NA"}' > ${matrix}_filler_${IntClust}

done


#create final big matrix
pastelist=$(for IntClust in {1..3} 4- 4+ {5..10} NA; do
    echo ${matrix}_ready_${IntClust}_CNA
    echo ${matrix}_filler_${IntClust}
  done)
#add gene list
cat $genes | awk 'BEGIN{print "\t"}{print}' > ${genes}_list

paste ${genes}_list $pastelist > ${matrix}_FINAL



#do separately, to determine gene order:
#sort the genes according to number of samples with amplification or deletion
cat ${matrix}_FINAL | \
  awk -vOFS="\t" '(NR>1){
    gene[NR]=$1
    for(i=2;i<=NF;i++)
      c[$1"-"$i]++
  }
  END{
    for(i=2;i<=NR;i++)
    print gene[i],c[gene[i]"-"0]+0,c[gene[i]"-"1]+0,c[gene[i]"-"2]+0,c[gene[i]"-NA"]+0,c[gene[i]"-"3],c[gene[i]"-"1]-c[gene[i]"-"2]+0
  }' | sort -k7,7nr | less
#use this sorting schema,
#then group together those that are close
ERBB2    412  105  16  24  27  89
CCND1    481  51   1   24  27  50
MYC      483  49   1   24  27  48
RPS6KB1  501  26   6   24  27  20
ZNF703   453  50   30  24  27  20
PTRH2    502  25   6   24  27  19
APPBP2   502  24   7   24  27  17
RSF1     491  29   13  24  27  16
INTS4    488  30   15  24  27  15
PPM1D    505  21   7   24  27  14
PIK3CA   527  4    2   24  27  2
TBL1XR1  527  4    2   24  27  2
EMSY     500  17   16  24  27  1
GATA3    527  1    5   24  27  -4
KRAS     517  2    14  24  27  -12
CDKN1B   512  3    18  24  27  -15
AKT1     515  0    18  24  27  -18
RB1      486  0    47  24  27  -47
CDH1     466  0    67  24  27  -67
TP53     466  0    67  24  27  -67
CBFB     462  0    71  24  27  -71
MAP2K4   459  0    74  24  27  -74
PPP2R2A  434  0    99  24  27  -99


echo $matrix
/tmp/tmp.icm6EN21f8
#paste this temporary file name into the variable holder in R







########
