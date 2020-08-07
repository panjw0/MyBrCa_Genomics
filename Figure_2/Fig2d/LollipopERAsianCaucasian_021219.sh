


#######################################################
#grab mutations to be displayed using lollipop plots
#compare ERpos/ERneg in Asians/Caucasians
#######################################################

maf=$(mktemp)
samples=$(mktemp)

####################################
#             MyBrCa
####################################
cd ~/work/tumor_project/data/redo_samples/9_Strelka2_indels/annotated_MAF
cat ../../../../annotations/Study_metadata_8.5.19.csv | \
  awk -vFS="," 'NR>1 && $8=="P"{print $1}' > ${samples}_MyBrCa_ERpos
cat ../../../../annotations/Study_metadata_8.5.19.csv | \
  awk -vFS="," 'NR>1 && $8=="N"{print $1}' > ${samples}_MyBrCa_ERneg

#create maf small files
for ERstat in ERpos ERneg; do
  cat compiled_clean.maf | \
    awk -vFS="\t" -vp=${samples}_MyBrCa_${ERstat} \
      'BEGIN{
        while((getline<p)>0) t[$1]=1
      }
      ($1!~/^#/) && ($1!="Hugo_Symbol") && ($16 in t)' > ${maf}_MyBrCa_${ERstat}_clean
done


####################################
#            METABRIC
####################################
cat /home/mzabidi/work/tumor_project/data/published_data/brca_metabric/vcfs/annotated_MAF/compiled_clean.maf | \
  awk -vp=$gene -vFS="\t" '$1!~"#" && $1!="Hugo_Symbol" && $1==p && $17~/^MB/' > ${maf}_METABRIC_all_clean
cat /home/mzabidi/work/tumor_project/data/published_data/brca_metabric/data_clinical_patient.txt | \
  awk -vFS="\t" '$1!~/^#/ && $1!="PATIENT_ID" && $7=="Negative"{print $1}' > ${samples}_METABRIC_ERneg
cat /home/mzabidi/work/tumor_project/data/published_data/brca_metabric/data_clinical_patient.txt | \
  awk -vFS="\t" '$1!~/^#/ && $1!="PATIENT_ID" && $7=="Positve"{print $1}' > ${samples}_METABRIC_ERpos

#create maf small files
for ERstat in ERpos ERneg; do
  cat ${maf}_METABRIC_all_clean | \
    awk -vFS="\t" -vp=${samples}_METABRIC_${ERstat} \
      'BEGIN{
        while((getline<p)>0) t[$1]=1
      }
      ($17 in t)' > ${maf}_METABRIC_${ERstat}_clean
done


####################################
#  NikZainal, Korean and NonKorean
####################################
cd ~/work/tumor_project/data/published_data/NikZainal_291018/nextera_regions
awk -vFS="," '$1~/PD2/ && $(NF-2)=="negative"{print $1}' ../SuppTable1.csv | \
  sort | uniq -c | sort -k1,1nr | awk '{print $2}' > ${samples}_NikZainal_Korean_ERneg
awk -vFS="," '$1~/PD2/ && $(NF-2)=="positive"{print $1}' ../SuppTable1.csv | \
  sort | uniq -c | sort -k1,1nr | awk '{print $2}' > ${samples}_NikZainal_Korean_ERpos
awk -vFS="," '$1!~/PD2/ && $(NF-2)=="negative"{print $1}' ../SuppTable1.csv | \
  sort | uniq -c | sort -k1,1nr | awk '{print $2}' > ${samples}_NikZainal_Nonkorean_ERneg
awk -vFS="," '$1!~/PD2/ && $(NF-2)=="positive"{print $1}' ../SuppTable1.csv | \
  sort | uniq -c | sort -k1,1nr | awk '{print $2}' > ${samples}_NikZainal_Nonkorean_ERpos

#create maf small files
for race in Korean Nonkorean; do
  for ERstat in ERpos ERneg; do
    cat annotated_MAF/compiled_clean.maf | \
      awk -vFS="\t" -vp=${samples}_NikZainal_${race}_${ERstat} \
        'BEGIN{
          while((getline<p)>0)
          {
            t[$1"a"]=1
            t[$1"a2"]=1
          }
        }
        ($1!~/^#/) && ($1!="Hugo_Symbol") && ($16 in t)' > ${maf}_NikZainal_${race}_${ERstat}_clean
  done
done

####################################
#  TCGA, Asian and Caucasians
####################################
cd ~/work/tumor_project/data/published_data/brca_tcga
grep -v "#" data_bcr_clinical_data_patient.txt | \
  awk -vFS="\t" '$9=="ASIAN" && $1!="OTHER_PATIENT_ID"{print $2}' | \
    awk -vp=TCGA_EXCLUDE.txt 'BEGIN{while((getline<p)>0) if(++l>1) k[$1]=1}
      !($1 in k)' > ${samples}_TCGA_Asian
grep -v "#" data_bcr_clinical_data_patient.txt | \
  awk -vFS="\t" '$9=="WHITE" && $1!="OTHER_PATIENT_ID"{print $2}' | \
    awk -vp=TCGA_EXCLUDE.txt 'BEGIN{while((getline<p)>0) if(++l>1) k[$1]=1}
      !($1 in k)' > ${samples}_TCGA_Caucasian

cat SuppTable1.csv | awk -vFS="," 'NR>1 && ($4=="Negative"){print $1}' | \
  sort | uniq > ${samples}_TCGA_ERneg
cat SuppTable1.csv | awk -vFS="," 'NR>1 && ($4=="Positive"){print $1}' | \
  sort | uniq > ${samples}_TCGA_ERpos

#create maf small files
for race in Asian Caucasian; do
  for ERstat in ERpos ERneg; do
    cat vcfs/annotated_MAF/compiled_clean.maf | \
      awk -vFS="\t" -vp=${samples}_TCGA_${race} \
        'BEGIN{
          while((getline<p)>0) t[$1"-01"]=1
        }
        ($1!~/^#/) && ($1!="Hugo_Symbol") && ($16 in t)' | \
        awk -vFS="\t" -vp=${samples}_TCGA_${ERstat} \
          'BEGIN{
            while((getline<p)>0) t[$1"-01"]=1
          }
          ($16 in t)' > ${maf}_TCGA_${race}_${ERstat}_clean
  done
done

####################################
#            ZhengyanKan
####################################
cd ~/work/tumor_project/data/published_data/ZhengyanKan2018
cat Supplementary_Data_1_Korean_only.csv | \
  awk -vFS="," 'NR>1 && ($7=="HER2+") || ($7=="TN"){print $1}' > ${samples}_ERneg
cat Supplementary_Data_1_Korean_only.csv | \
  awk -vFS="," 'NR>1 && ($7=="ER+HER2+") || ($7=="ER+"){print $1}' > ${samples}_ERpos
#create maf small files
for ERstat in ERpos ERneg; do
  cat annotated_MAF/compiled_clean.maf | \
    awk -vOFS="\t" -vFS="\t" -vp=${samples}_${ERstat} \
      'BEGIN{
        while((getline<p)>0) t[$1]=1
      }
      $1!~"#" && $1!~"Hugo_Symbol" && ($16 in t)' > ${maf}_ZhengyanKan_${ERstat}_clean
done


##############################################
#   Prepare file for Mutation Mapper input
##############################################
gene=TP53
#also do for GATA3
gene=GATA3
for ERstat in ERpos ERneg; do
  #Asians
  cat ${maf}_MyBrCa_${ERstat}_clean ${maf}_ZhengyanKan_${ERstat}_clean \
    ${maf}_TCGA_Asian_${ERstat}_clean ${maf}_NikZainal_Korean_${ERstat}_clean | \
    awk -vp=$gene '$1==p' | \
      awk -vFS="\t" -vOFS="\t" 'BEGIN{
        l1="Sample_ID\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position"
        l2="Reference_Allele\tVariant_Allele\tMutation_Type\tProtein_Change"
        print l1,l2}
        {print $16,$1,$5,$6,$7,$12,$13,$9,$42}' | awk '!seen[$0]++' > ~/R_output/Asians_${gene}_${ERstat}.txt
  #Caucasians
  cat ${maf}_METABRIC_${ERstat}_clean ${maf}_TCGA_Caucasian_${ERstat}_clean ${maf}_NikZainal_Nonkorean_${ERstat}_clean | \
    awk -vp=$gene '$1==p' | \
      awk -vFS="\t" -vOFS="\t" 'BEGIN{
        l1="Sample_ID\tHugo_Symbol\tChromosome\tStart_Position\tEnd_Position"
        l2="Reference_Allele\tVariant_Allele\tMutation_Type\tProtein_Change"
        print l1,l2}
        {print $16,$1,$5,$6,$7,$12,$13,$9,$42}' | awk '!seen[$0]++' > ~/R_output/Caucasians_${gene}_${ERstat}.txt

done

for ERstat in ERpos ERneg; do
  echo $ERstat
  for race in Asians Caucasians; do
    echo $race $(awk 'NR>1' ~/R_output/${race}_${gene}_${ERstat}.txt | wc -l) \
      $(awk 'NR>1{print $1}' ~/R_output/${race}_${gene}_${ERstat}.txt | sort | uniq | wc -l) \
      $(awk 'NR>1{print $9}' ~/R_output/${race}_${gene}_${ERstat}.txt | sort | uniq | wc -l)
  done
done
#TP53
#total_records, samples, position-aachanges
ERpos
Asians 165 164 117
Caucasians 467 454 255
ERneg
Asians 224 223 119
Caucasians 528 518 260




#use annotations from Mutation Mapper,
#that refers to PFAM annotation
#count #mutations in transactivating (a.a 6-29; chrpos 7661779-7687550),
#DNA-binding (95-288),
#and tetramerization(318-358) domains

#found no good cDNA->a.a simple conversion
#do it manually
#TP53 runs in the reverse direction

#Caucasian ERneg
ERstat=ERneg
#transactivating: 2 truncating
#DNA-binding:
#from MB-7159/T102fs till MB-0613/E287*
sort -k4,4nr ~/R_output/Caucasians_${gene}_${ERstat}.txt | \
  awk '$1!="Sample_ID"' | \
    awk '($1=="MB-7159" && $9=="p.T102fs"){y=1}(y)' | \
    awk '{
      if(!($1=="MB-0613" && $9=="p.E287*")){print}
    else{print ; exit}
    }' | awk '{print $8}' | sort | uniq -c
    54 Frame_Shift_Del
    16 Frame_Shift_Ins
    12 In_Frame_Del
     1 In_Frame_Ins
   270 Missense_Mutation
    53 Nonsense_Mutation
    32 Splice_Site
#tetramerization:
#from MB-6246/P322fs till MB-4757/E349* and MB-7020/ELKDAQAG349del
#till the end basically, so it's OK
sort -k4,4nr ~/R_output/Caucasians_${gene}_${ERstat}.txt | \
  awk '$1!="Sample_ID"' | \
    awk '($1=="MB-6246" && $9=="p.P322fs"){y=1}(y)' | \
      awk '{print $8}' | sort | uniq -c
      9 Frame_Shift_Del
      1 Frame_Shift_Ins
      1 In_Frame_Del
     11 Missense_Mutation
     13 Nonsense_Mutation
      5 Splice_Site

#Caucasian ERpos
ERstat=ERpos
#transactivating: 2 missense
#DNA-binding:
#from PD7238a/ till TCGA-BH-A0DL-01/E286A
sort -k4,4nr ~/R_output/Caucasians_${gene}_${ERstat}.txt | \
  awk '$1!="Sample_ID"' | \
    awk '($1=="PD7238a" && $9=="p.SV96fs"){y=1}(y)' | \
    awk '{
      if(!($1=="TCGA-BH-A0DL-01" && $9=="p.E286A")){print}
    else{print ; exit}
    }' | awk '{print $8}' | sort | uniq -c
    35 Frame_Shift_Del
    13 Frame_Shift_Ins
    11 In_Frame_Del
     3 In_Frame_Ins
   273 Missense_Mutation
    38 Nonsense_Mutation
    23 Splice_Site
#tetramerization:
#from PD9760a/K321* till MB-7263/D352Y
sort -k4,4nr ~/R_output/Caucasians_${gene}_${ERstat}.txt | \
  awk '$1!="Sample_ID"' | \
    awk '($1=="PD9760a" && $9=="p.K321*"){y=1}(y)' | \
    awk '{
      if(!($1=="MB-7263" && $9=="p.D352Y")){print}
    else{print ; exit}
    }' | awk '{print $8}' | sort | uniq -c
5 Frame_Shift_Del
2 Frame_Shift_Ins
1 In_Frame_Del
12 Missense_Mutation
11 Nonsense_Mutation
4 Splice_Site





#Asian ERneg
ERstat=ERneg
#transactivating: NO mutations
#DNA-binding:
#from SD0986/P98fs till SD1530/E286Q
sort -k4,4nr ~/R_output/Asians_${gene}_${ERstat}.txt | \
  awk '$1!="Sample_ID"' | \
    awk '($1=="SD0986" && $9=="p.P98fs"){y=1}(y)' | \
    awk '{
      if(!($1=="SD1530" && $9=="p.E286Q")){print}
    else{print ; exit}
    }' | awk '{print $8}' | sort | uniq -c
    18 Frame_Shift_Del
     7 Frame_Shift_Ins
     4 In_Frame_Del
   137 Missense_Mutation
    17 Nonsense_Mutation
    11 Splice_Site
#tetramerization:
#from PD22363a-splice till TCGA-C8-A12Z-01/p.R342P
#till the end basically, so it's OK
sort -k4,4nr ~/R_output/Asians_${gene}_${ERstat}.txt | \
  awk '$1!="Sample_ID"' | \
    awk '($1=="PD22363a" && $8=="Splice_Site"){y=1}(y)' | \
      awk '{print $8}' | sort | uniq -c
    1 Frame_Shift_Del
    3 Missense_Mutation
    4 Nonsense_Mutation
    3 Splice_Site






#Asian ERpos
ERstat=ERpos
#transactivating: 1 truncating
#DNA-binding:
#from SD0891/Q100fs till SD0164/p.E285*
sort -k4,4nr ~/R_output/Asians_${gene}_${ERstat}.txt | \
  awk '$1!="Sample_ID"' | \
    awk '($1=="SD0891" && $9=="p.Q100fs"){y=1}(y)' | \
    awk '{
      if(!($1=="SD0164" && $9=="p.E285*")){print}
    else{print ; exit}
    }' | awk '{print $8}' | sort | uniq -c
    13 Frame_Shift_Del
     3 Frame_Shift_Ins
     3 In_Frame_Del
    96 Missense_Mutation
    15 Nonsense_Mutation
     7 Splice_Site
#tetramerization:
#from SD0916/L323fs till BR408p.R342*
#(till except the very last one)
sort -k4,4nr ~/R_output/Asians_${gene}_${ERstat}.txt | \
  awk '$1!="Sample_ID"' | \
    awk '($1=="SD0916" && $9=="p.L323fs"){y=1}(y)' | \
    awk '{
      if(!($1=="BR408" && $9=="p.R342*")){print}
    else{print ; exit}
    }' | awk '{print $8}' | sort | uniq -c
    3 Frame_Shift_Del
    1 Missense_Mutation
    5 Nonsense_Mutation
    1 Splice_Site










#########
