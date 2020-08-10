


##########################################################################
#   heatmap of mutational signature weights vs immune score for TCGA
##########################################################################
cd ~/work/tumor_project/data/published_data/brca_tcga/vcfs

#grab immune scores and append to signature weights
for race in Caucasian Asian; do
  cat signatures_table_${race}s.txt | \
    awk -vOFS="\t" -vp=/home/mzabidi/work/tumor_project/annotations/immune_score/TCGA_immunescores.txt \
    'BEGIN{
      while((getline<p)>0)
      {
        IMPRES[$1]=$5
        eIFNg[$1]=$6
        GSVABindea[$1]=$7
        ESTIMATE[$1]=$8
      }
    }
    (NR==1){print $0,"IMPRES","eIFNgamma","GSVABindea","ESTIMATE"}
    (NR>1){print $0,IMPRES[$1]?IMPRES[$1]:"NA",eIFNg[$1]?eIFNg[$1]:"NA",GSVABindea[$1]?GSVABindea[$1]:"NA",ESTIMATE[$1]?ESTIMATE[$1]:"NA"}' | \
     grep -v NA > ${x}_grandtable_immune_${race}
done

echo $x
/tmp/tmp.z5DFKRsdXU






########
