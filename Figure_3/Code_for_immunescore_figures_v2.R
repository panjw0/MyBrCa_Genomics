### Code for Immune Score Figures ###


#### Fig 4A ###

### ESTIMATE score ###

df2 <- read.table("ESTIMATE_score_R_input_2.txt", header=TRUE)
df3<-df2[df2$Cohort!="TCGA_Asn",]
df3$Cohort<-factor(df3$Cohort, levels=c("MyBrCa","Kan","TCGA_Cau"))



par(mar=c(4,6,2,2))
boxplot(ImmuneScore~Cohort, data=df3, col=c("Plum","Purple","Red"), outline=F, frame=F,ylab=NULL)
title(ylab="ESTIMATE Immune Score", line=4, cex.lab=1.2)

fit2 <- aov(ImmuneScore~Cohort, data=df3)
summary(fit2)
TukeyHSD(fit2)




### Bindea score ###


df2 <- read.table("Bindea_score_R_input_2.txt", header=TRUE)
df3<-df2[df2$Cohort!="TCGA_Asn",]
df3$Cohort<-factor(df3$Cohort, levels=c("MyBrCa","Kan","TCGA_Cau"))

par(mar=c(4,6,2,2))
boxplot(Bindea_GSVA~Cohort, data=df3, col=c("Plum","Purple","Red"), outline=F, frame=F, ylab=NULL)
title(ylab="Bindea Immune Score", line=4, cex.lab=1.2)

fit2 <- aov(Bindea_GSVA~Cohort, data=df3)
summary(fit2)
TukeyHSD(fit2)


### EIFNg score ###


df2 <- read.table("EIFNg_score_R_input_2.txt", header=TRUE)
df3<-df2[df2$Cohort!="TCGA_Asn",]
df3$Cohort<-factor(df3$Cohort, levels=c("MyBrCa","Kan","TCGA_Cau"))

par(mar=c(4,6,2,2))
boxplot(Expanded_IFNg~Cohort, data=df3, col=c("Plum","Purple","Red"), outline=F, frame=F,ylab=NULL)
title(ylab="Expanded IFN-gamma Score", line=4, cex.lab=1.2)

fit2 <- aov(Expanded_IFNg~Cohort, data=df3)
summary(fit2)
TukeyHSD(fit2)


## IMPRES Score ##

df2 <- read.table("IMPRES_score_R_input_2.txt", header=TRUE)
df3<-df2[df2$Cohort!="TCGA_Asn",]
df3$Cohort<-factor(df3$Cohort, levels=c("MyBrCa","Kan","MB","TCGA_Cau"))

par(mar=c(4,4,2,2))
boxplot(IMPRES_score~Cohort, data=df3, col=c("Plum","Purple","Brown","Red"), outline=F, frame=F)
title(ylab="ESTIMATE Immune Score", line=4, cex.lab=1.2)

fit2 <- aov(IMPRES_score~Cohort, data=df3)
summary(fit2)
TukeyHSD(fit2)



### CIBERSORT CD8 boxplot ###


df2 <- read.table("CIBERSORT_R_input_2.txt", header=TRUE)
df3<-df2[df2$Cohort!="TCGA_Asn",]
df3$Cohort<-factor(df3$Cohort, levels=c("MyBrCa","Kan","TCGA_Cau"))

par(mar=c(4,6,2,2))
boxplot(T.cells.CD8~Cohort, data=df3, col=c("Plum","Purple","Red"), outline=F, frame=F)#, las=2)





##### Figure 4B IMPRES vs Cohort across subtypes ###



df2 <- read.table("IMPRES_v_cohort_across_ic10_subtypes_R_input_2.txt", header=TRUE)


df3<-df2[df2$Cohort=="TCGA_Cau"|df2$Cohort=="MyBrCa",]
df3$IntClust2_Cohort<-factor(df3$IntClust2_Cohort, levels=c("1_MyBrCa","1_TCGA",
                                                            "2_MyBrCa","2_TCGA",
                                                            "3_MyBrCa","3_TCGA",
                                                            "4+_MyBrCa","4+_TCGA",
                                                            "4-_MyBrCa","4-_TCGA",
                                                            "5_MyBrCa","5_TCGA",
                                                            "6_MyBrCa","6_TCGA",
                                                            "7_MyBrCa","7_TCGA",
                                                            "8_MyBrCa","8_TCGA",
                                                            "9_MyBrCa","9_TCGA",
                                                            "10_MyBrCa","10_TCGA"))



par(mar=c(8,8,2,2))
boxplot(IMPRES_score~IntClust2_Cohort, data=df3, col=c("Plum","Red"), outline=F, frame=F, las=2, ylab=NULL, xlab=NULL)
title(ylab="IMPRES Score", line=3, cex.lab=1.2)



t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="1",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="2",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="3",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="4+",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="4-",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="5",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="6",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="7",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="8",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="9",])
t.test(IMPRES_score~Cohort, data=df3[df3$IntClust2=="10",])




###### IHC Validation scatterplot #####


df2 <- read.table("ImmuneScore_IHC_validation.txt", header=TRUE)

plot(IMPRES_score~CD8_pct, data=df2, 
     xlab="CD8 Staining (% area)", ylab="IMPRES Score", pch=1)


abline(lm(IMPRES_score~CD8_pct, data=df2), col="red") 






