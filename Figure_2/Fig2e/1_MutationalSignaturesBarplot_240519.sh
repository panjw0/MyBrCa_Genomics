


##########################################################
#draw barplots of %samples carrying signature X
##########################################################
cd ~/work/tumor_project/data

#Asians
SJMC=/people/mzabidi/work/tumor_project/data/redo_samples/9_Strelka2_indels/deconstructSigs_signatures_breast_only/signatures_table.txt
SMC=/people/mzabidi/work/tumor_project/data/published_data/ZhengyanKan2018/vcfs/deconstructSigs_signatures/signatures_table.txt
TCGA_Asian=/people/mzabidi/work/tumor_project/data/published_data/brca_tcga/vcfs/signatures_table_Asians.txt
NikZainal_Korean=/people/mzabidi/work/tumor_project/data/published_data/NikZainal_291018/nextera_regions/deconstructSigs_signatures_breast_sigs/signatures_table_Koreans.txt

#Caucasian
TCGA_Caucasian=/people/mzabidi/work/tumor_project/data/published_data/brca_tcga/vcfs/signatures_table_Caucasians.txt
NikZainal_NonKorean=/people/mzabidi/work/tumor_project/data/published_data/NikZainal_291018/nextera_regions/deconstructSigs_signatures_breast_sigs/signatures_table_NonKoreans.txt

cutoff=0.2
datasets=( $SJMC $SMC $TCGA_Asian $NikZainal_Korean $TCGA_Caucasian $NikZainal_NonKorean )
suffix=( sjmc smc tcga_asian nikzainal_korean tcga_caucasian nikzainal_nonkorean )

table=$(mktemp)
for i in ${!datasets[*]}; do
  cat ${datasets[$i]} | \
    awk -vOFS="\t" -vp=$cutoff \
      'NR==1{
        for(i=1;i<=NF;i++)
          sig[i]=$i
        tot_all_sigs=NF
      }
      NR>1{
        for(i=1;i<=NF;i++)
          if($(i+1)>=p) count_sig[i+1]++
      }
      END{
        for(i=1;i<=tot_all_sigs;i++)
          print sig[i],count_sig[i+1]+0,NR-1
      }' > ${table}_${suffix[$i]}
done

echo $table
/tmp/tmp.CSssKCHzdS






###########
