library(survival)
library(survminer)
library(dplyr)
library("RColorBrewer")



cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


sa.df<-read.table('MyBrCa_survival_data.txt',header=T)
sa.df$time_int_yrs<-sa.df$time_int/365.25
df<-sa.df


### Plot unadjusted with log rank test ###



## Figure 5a: ER status


sa.df<-df
surv_object <- Surv(time = sa.df$time_int_yrs, event = sa.df$Death1)
fit1 <- survfit(surv_object ~ ER, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)



## Figure 5a: IntClust



sa.df<-df[df$IntClust_11=="3"|df$IntClust_11=="5"|df$IntClust_11=="8"|df$IntClust_11=="10",]
surv_object <- Surv(time = sa.df$time_int_yrs, event = sa.df$Death1)
fit1 <- survfit(surv_object ~ IntClust_11, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)


## Figure 5a: PAM50



sa.df<-df[df$PAM50!="Normal",]
surv_object <- Surv(time = sa.df$time_int_yrs, event = sa.df$Death1)
fit1 <- survfit(surv_object ~ PAM50, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)



## Figure 5b: All patients



sa.df<-df
surv_object <- Surv(time = sa.df$time_int_yrs, event = sa.df$Death1)
fit1 <- survfit(surv_object ~ TP53_sm, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)
fit1 <- survfit(surv_object ~ IMPRES_group, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)
fit1 <- survfit(surv_object ~ TP53_sm+IMPRES_group, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)



## Figure 5b: ER+ patients



sa.df<-df[df$ER=="P",]
surv_object <- Surv(time = sa.df$time_int_yrs, event = sa.df$Death1)
fit1 <- survfit(surv_object ~ TP53_sm, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)
fit1 <- survfit(surv_object ~ IMPRES_group, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)
fit1 <- survfit(surv_object ~ TP53_sm+IMPRES_group, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)



## Figure 5b: ER- patients


sa.df<-df[df$ER=="N",]
surv_object <- Surv(time = sa.df$time_int_yrs, event = sa.df$Death1)
fit1 <- survfit(surv_object ~ TP53_sm, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)
fit1 <- survfit(surv_object ~ IMPRES_group, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)
fit1 <- survfit(surv_object ~ TP53_sm+IMPRES_group, data = sa.df)
ggsurvplot(fit1, data = sa.df, pval = TRUE, xlim=c(0,5),palette=cbbPalette)





#### Plot adjusted cox proportional hazard model ###


sa.df<-df
cfit4a <- coxph(Surv(time = sa.df$time_int_yrs, event = sa.df$Death1) ~ GradeIndex+AgeDiag1+Stage_group+ER+IMPRES_group+TP53_sm,
                data=sa.df)


ggforest(cfit4a, data = sa.df)

