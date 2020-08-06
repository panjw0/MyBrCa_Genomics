
#### Code for Figure 1: Subtypes ####


df2 <- read.table("Subtypes_R_input_7.txt", header=TRUE)
df2$Cohort<-factor(df2$Cohort, levels=c("MyBrCa","Kan","TCGA_Asn","MB","TCGA_Cau"))
df2$Cohort2<-factor(df2$Cohort2, levels=c("MyBrCa","Asian","Caucasian"))
df2$PAM50<-factor(df2$PAM50, levels=c("LumA","LumB","Her2","Basal","Normal"))
df2$IntClust2<-factor(df2$IntClust2, levels=c("1","2","3","4+","4-","5","6","7","8","9","10"))


# Grouped Bar Plot




cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




#### Fig 1A ###

counts <- table(df2$Cohort2, df2$IntClust2)
counts2 <- t(counts)
counts3 <- scale(counts2, FALSE, colSums(counts2)) * 100
counts4<-t(counts3)



barplot(counts4,
        xlab="IntClust Subtype", col=cbPalette[1:3],
        ylab="Frequency (%)",
        beside=TRUE)

legend("right", 
       legend = c("MyBrCa","Other Asian","Caucasian"), 
       fill = cbPalette[1:3], ncol = 1,
       cex = 0.7)






### Fig 1B ###



df3<-df2[df2$Menopause_inferred_50=="Pre",]

counts <- table(df3$Cohort2, df3$IntClust2)
counts2 <- t(counts)
counts3 <- scale(counts2, FALSE, colSums(counts2)) * 100
counts4<-t(counts3)




barplot(counts4,
        xlab="IntClust Subtype for Premenopausal Women (Age <50)", col=cbPalette[1:3],
        ylab="Frequency (%)",
        beside=TRUE)

legend("right", 
       legend = c("MyBrCa","Other Asian","Caucasian"), 
       fill = cbPalette[1:3], ncol = 1,
       cex = 0.7)






### Fig 1C ###

df3<-df2[df2$Menopause_inferred_50=="Post",]

counts <- table(df3$Cohort2, df3$IntClust2)
counts2 <- t(counts)
counts3 <- scale(counts2, FALSE, colSums(counts2)) * 100
counts4<-t(counts3)




barplot(counts4,
        xlab="IntClust Subtype for Postmenopausal Women (Age >50)", col=cbPalette[1:3],
        ylab="Frequency (%)",
        beside=TRUE)

legend("right", 
       legend = c("MyBrCa","Other Asian","Caucasian"), 
       fill = cbPalette[1:3], ncol = 1,
       cex = 0.7)






### Chi-square tests ###


#example: IC5

df4<-df2
counts <- table(df4$Race, df4$IntClust2=="5")
counts.m<-as.matrix(counts)
chisq.test(counts.m)






