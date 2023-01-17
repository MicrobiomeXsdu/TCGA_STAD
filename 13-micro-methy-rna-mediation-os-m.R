#BiocManager::install("timereg")
library(glmnet)
library(caret)
library(timereg)
library(survival)
data(ses_dat_m)
###############加载数据并整理，最终为data_all,micro_fil和methy_fil###########
setwd("/home/yuekaile/keller/stad_rna/micro-methy-bind-os//")
micro <- read.csv("/home/yuekaile/keller/stad_rna/micro_os/coxph-stad-os-micro-data-13.csv",
                  header = T,row.names = 1,check.names = F)
methy_dmp <- read.csv("/home/yuekaile/keller/ykl/methy-os/coxph-stad-os-methy-dmp-data-12.csv",
                      header = T,row.names = 1,check.names = F)
methy_dmr <- read.csv("/home/yuekaile/keller/ykl/methy-os/coxph-stad-os-methy-dmr-data-12.csv",
                      header = T,row.names = 1,check.names = F)
meta <- read.csv("/home/yuekaile/keller/stad_rna/pheno_fil.csv",
               header = T,row.names = 1,check.names = F)
meta2 <- read.csv("/home/yuekaile/keller/tcga/TCGA-STAD-methylation-phenotype.csv",
                 header = T,row.names = 1,check.names = F)

micro$age <- meta2$age_at_initial_pathologic_diagnosis[match(rownames(micro),rownames(meta2))]
micro$gender <- meta2$gender.demographic[match(rownames(micro),rownames(meta2))]
micro$ostime <- meta$ostime[match(rownames(micro),meta$submitter_id_samples)]

micro_fil <- micro[rownames(methy_dmp),]
all(which(rownames(micro_fil)==rownames(methy_dmp)))
all(which(rownames(micro_fil)==rownames(methy_dmr)))
which(colnames(methy_dmp)==colnames(methy_dmr))

methy_fil <- cbind(methy_dmp,methy_dmr)
 
data_all <- cbind(methy_fil,micro_fil)
agea <- data_all$age
genderr <- data_all$gender

#拟合模型
for (i in 1:57) {
  x <- data_all[,1]
  for (j in 58:70) {
    m <- data_all[,65]
  ols_m <- glm(m~ x + agea + factor(genderr),
             data=data_all)
  aa <- summary(ols_m)
  print(c(i,j,aa[["coefficients"]][2,1],aa[["coefficients"]][2,4]))
  aalen_m2 <- aalen(Surv(ostime,os) ~ x 
                  + m
                  + agea
                  + const(factor(genderr)), 
                  data = data_all, robust=T ,n.sim = 100)
  obj <- try(p1 <- data.frame(aalen_m2[["pval.testBeqC"]][2],aalen_m2[["pval.testBeqC"]][3],aalen_m2[["pval.testBeqC"]][6]),silent=TRUE)
  if (is(obj, "try-error")){
    next
  } else if(aalen_m2[["pval.testBeqC"]][2] <0.05 & aalen_m2[["pval.testBeqC"]][3] <0.05 & aalen_m2[["pval.testBeqC"]][6] <0.05){
   print(summary(aalen_m2))
    
   }else{
    next
  }
  }
  }

for (i in 1:57) {
  m <- data_all[,i]
  for (j in 58:70) {
    x <- data_all[,j]
    ols_m <- glm(m~ x + agea + factor(genderr),
                 data=data_all)
    aa <- summary(ols_m)
    print(c(j,i,aa[["coefficients"]][2,1],aa[["coefficients"]][2,4]))
    aalen_m2 <- aalen(Surv(ostime,os) ~ x 
                      + m+m*x
                      + agea
                      + const(factor(genderr)), 
                      data = data_all, robust=T ,n.sim = 100)
    obj <- try(p1 <- data.frame(aalen_m2[["pval.testBeqC"]][2],aalen_m2[["pval.testBeqC"]][3],aalen_m2[["pval.testBeqC"]][6]),silent=TRUE)
    if (is(obj, "try-error")){
      next
    } 
    else if(aalen_m2[["pval.testBeqC"]][2] <0.05 & aalen_m2[["pval.testBeqC"]][3] <0.05 & aalen_m2[["pval.testBeqC"]][6] <0.05){
        print(summary(aalen_m2))
      }else{
        next
      }
    }
  }


# Additive Aalen Model 
# 
# Test for nonparametric terms 
# 
# Test for non-significant effects 
# Supremum-test of significance p-value H_0: B(t)=0
# (Intercept)                          3.78                0.02
# x                                    3.82                0.02
# m                                    2.45                0.19
# 
# Test for time invariant effects 
# Kolmogorov-Smirnov test p-value H_0:constant effect
# (Intercept)                         41.70                        0.07
# x                                   44.30                        0.07
# m                                    0.95                        0.06
# Cramer von Mises test p-value H_0:constant effect
# (Intercept)                       1320000                        0.06
# x                                 1460000                        0.06
# m                                     557                        0.07
# 
# Parametric terms :     
#   Coef.       SE Robust SE    z  P-val lower2.5% upper97.5%
#   const(agea)                1.32e-05 5.22e-06  5.48e-06 2.41 0.0160  2.97e-06   2.34e-05
# const(factor(genderr))male 2.30e-04 1.15e-04  1.13e-04 2.04 0.0418  4.60e-06   4.55e-04
# 
# Call: 
#   aalen(formula = Surv(ostime, os) ~ x + m + const(agea) + const(factor(genderr)), 
#         data = data_all, robust = T, n.sim = 100)
# 
# NULL



#################生存曲线##############
library(survminer) # 加载包
library(survival)
Surv(data_all$ostime,data_all$os)
data_all$risk <- as.factor(data_all$risk)
fit <- survfit(Surv(ostime,os) ~ risk,  # 创建生存对象 
               data = data_all) # 数据集来源
fit # 查看拟合曲线信息

summary(fit)

ggsurvplot(fit, data = lung)

ggsurvplot(fit, # 创建的拟合对象
           data = data_all,  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           xlab = "Follow up time(d)", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "", # 设置图例标题
           legend.labs = c("high", "low"), # 指定图例分组标签
           break.x.by = 300)  # 设置x轴刻度间距
