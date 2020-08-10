## Pyclone vs ESTIMATE groups scatterplot ##

df2 <- read.table("log_pyclone_vs_estimate.txt", header=TRUE)
df2<-df2[order(df2$ESTIMATE_ImmuneScore),]

group<-seq(0,494,by=26)

df3<-NULL

for (i in 1:20) {

  x<-1+group[i]
  y<-26+group[i]
  df3$Median_Log_Pyclone[i] <- median(df2$Log_PyClone_Clusters[x:y]) 
  df3$Median_ESTIMATE_Score[i] <- median(df2$ESTIMATE_ImmuneScore[x:y])
  
}

df3<-as.data.frame(df3)



plot(Median_ESTIMATE_Score~Median_Log_Pyclone, data=df3, 
     xlab="Median log2(Pyclone clusters)", ylab="Median ESTIMATE Immune Score", pch=1)


abline(lm(Median_ESTIMATE_Score~Median_Log_Pyclone, data=df3), col="red") 


