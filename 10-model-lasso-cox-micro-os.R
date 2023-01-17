load("~/keller/stad_rna/phy_genus_species.Rdata")
phy_species_sub <- prune_samples(phy_species@sam_data[["sample_type_samples"]]=="Primary Tumor",phy_species)

phy_species_sub@sam_data[["os"]] <- as.factor(phy_species_sub@sam_data[["os"]])
phy_species_sub <- prune_samples(phy_species_sub@sam_data[["os"]]=="1"|phy_species_sub@sam_data[["os"]]=="0",phy_species_sub)
table(phy_species_sub@sam_data[["os"]])
# 0   1 
# 204 146 
library(pacman)
library(magrittr)
library(reshape2)
library(MicrobiotaProcess)
library(phyloseq)
alpha_species <- get_alphaindex(phy_species_sub)
alpha_species_plot <- ggbox(alpha_species, geom="violin",factorNames="os")+
  scale_fill_manual(values=c("#66C2A5","#F98D63"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
alpha_species_plot

alpha_species_data <- cbind(alpha_species@alpha,alpha_species@sampleda)
write.csv(alpha_species_data,file = "alpha_species_data_os.csv")

#######beta多样性+adonis\anosim##############################
pacman::p_load(tidyverse,ggrepel,vegan,ape,ggsignif,patchwork,multcomp)
otu <-phy_species_sub@otu_table
groups <- phy_species_sub@sam_data %>%
  as.list()
#bray距离
pcoa <- vegdist(decostand(otu, "norm"),method = "bray") %>% 
  pcoa(correction = "none", rn = NULL)
#unifrac距离
#library(picante)
#pcoa <- picante::unifrac(otu,tree) %>%
#  pcoa(correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
pcoadata <- data.frame(rownames(pcoa$vectors),
                       PC1,PC2,groups$os)
colnames(pcoadata) <-c("sample","PC1","PC2","group")
pcoadata$pc1_abs <- abs(pcoadata$PC1)
pcoadata$pc2_abs <- abs(pcoadata$PC2)
pcoadata$pc_sum <- pcoadata$pc1_abs+pcoadata$pc2_abs
write.csv(pcoadata,"pcoadata_os-species.csv")
pcoadata$group <- factor(pcoadata$group,levels=c("0","1"))
str(pcoadata)
yf <- pcoadata
yd1 <- yf %>% group_by(group) %>% summarise(Max = max(PC1))
yd2 <- yf %>% group_by(group) %>% summarise(Max = max(PC2))
yd1$Max <- yd1$Max + max(yd1$Max)*0.1
yd2$Max <- yd2$Max + max(yd2$Max)*0.1

res1 <- aov(PC1~group,data = pcoadata) %>% 
  glht(linfct=mcp(group="Tukey")) %>% cld(alpah=0.05)
res2 <- aov(PC2~group,data = pcoadata) %>% 
  glht(linfct=mcp(group="Tukey")) %>% cld(alpah=0.05)

test <- data.frame(PC1 = res1$mcletters$Letters,PC2 = res2$mcletters$Letters,
                   yd1 = yd1$Max,yd2 = yd2$Max,group = yd1$group)
test$group <- as.factor(test$group)


p1 <- ggplot(pcoadata, aes(PC1, PC2)) +
  geom_point(aes(colour=group),size=2.5)+
  scale_color_manual(values=c("#66C2A5","#F98D63"))+
  stat_ellipse(aes(colour=group),level=0.95,show.legend=F,linetype="dotted")+#aes(colour=group),
  scale_fill_manual(values=c("#66C2A5","#F98D63"))+
  labs(x=(floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PC1 ( ", ., "%", " )"),
       y=(floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PC2 ( ", ., "%", " )")) +
  theme(text=element_text(size=10))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        axis.title.x=element_text(colour='black', size=10),
        axis.title.y=element_text(colour='black', size=10),
        axis.text=element_text(colour='black',size=12),
        legend.title=element_blank(),
        legend.key.height=unit(0.6,"cm"),
        legend.position = c(0.75, 0.95),legend.direction = "horizontal")
p1
p2 <- ggplot(pcoadata,aes(group,PC1)) +
  geom_boxplot(aes(fill = group))+
  scale_fill_manual(values=c("#66C2A5","#F98D63"))+
  geom_jitter(shape=16,size=1.5,position=position_jitter(0.2))+
  geom_text(data = test,aes(x = group,y = yd1,label = PC1),
            size = 5,color = "black",fontface = "plain")+
  theme(panel.background = element_rect(fill='white',
                                        colour='black'))+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(colour='black',size=10,face = "plain"),
        axis.text.x=element_blank(),
        legend.position = "none")+coord_flip()
p2 

p3 <- ggplot(pcoadata,aes(group,PC2)) +
  geom_boxplot(aes(fill = group),levels=c("no","yes")) +
  scale_fill_manual(values=c("#66C2A5","#F98D63"))+
  geom_jitter(shape=16,size=1.5,position=position_jitter(0.2))+
  # geom_text(data = test,aes(x = group,y = yd2,label = PC2),
  #           size = 5,color = "black",fontface = "plain")+
  theme(panel.background = element_rect(fill='white',
                                        colour='black'))+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour='black',size=10,angle = 0,
                                 vjust = 1,hjust = 0.5,face = "plain"),
        axis.text.y=element_blank(),
        legend.position = "none")

p3

otu.adonis=adonis(otu~group,data = pcoadata,distance = "bray")
p4 <- ggplot(pcoadata,
             aes(PC1, PC2))+
  geom_text(aes(x = -0.5,
                y = 0.6,
                label = paste("PERMANOVA:\ndf = ",
                              otu.adonis$aov.tab$Df[1],"\nR2 = ",
                              round(otu.adonis$aov.tab$R2[1],4),
                              "\np-value = ",
                              otu.adonis$aov.tab$`Pr(>F)`[1],
                              sep = "")),size = 2.8) +theme_bw() +
  xlab(NULL) + ylab(NULL) +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
p4
p2+p4+p1+p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)        

####################anosim#########################
#https://blog.csdn.net/weixin_45822007/article/details/121623750?ops_request_misc=&request_id=&biz_id=102&utm_term=anosim代码&utm_medium=distribute.pc_search_result.none-task-blog-2~blog~sobaiduweb~default-3-121623750.nonecase&spm=1018.2226.3001.4450
set.seed(5555)
anosim_species <-anosim(phy_species_sub@otu_table,phy_species_sub@sam_data[["os"]],
                        permutations = 999, distance = "bray")
summary(anosim_species)
# ANOSIM statistic R: 0.02483 
# Significance: 0.002 

dc_species = data.frame(dis=anosim_species$dis.rank, class=anosim_species$class.vec)
# 多重比较
pairwise.anosim <-function(x,factors, sim.method, p.adjust.m){
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  R = c()
  p.value = c()  
  for(elem in 1:ncol(co)){
    ad = anosim(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),],                
                factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))], permutations = 999,distance = "bray");    
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));    
    R = c(R,ad$statistic);   
    p.value = c(p.value,ad$signif)  
  }
  p.adjusted =p.adjust(p.value,method=p.adjust.m) 
  pairw.res = data.frame(pairs, R, p.value,p.adjusted) 
  return(pairw.res)
} # 不需要更改
set.seed(5555)
pairwise.anosim(phy_species_sub@otu_table,phy_species_sub@sam_data[["os"]], sim.method="bray", p.adjust.m= "fdr")
# pairs          R p.value p.adjusted
# 1 no vs yes 0.06108336   0.026      0.026
set.seed(555)
pairwise.anosim(phy_species_sub@otu_table,phy_species_sub@sam_data[["pathologic_M"]], sim.method="bray", p.adjust.m= "fdr")
# pairs           R p.value p.adjusted
# 1 M0 vs MX  0.04777519   0.096      0.144
# 2 M0 vs M1  0.06473663   0.027      0.081
# 3 MX vs M1 -0.01857059   0.684      0.684

# * A=0.1881，A>0表示组间差异大于组内差异，A<0表示组间差异小于组内差异。
# * observed delta=0.2667,值越小说明组内差异越小。
# * expected delta=0.3285,值越小说明组间差异越小。
# * significance of delta(P)=0.001，小于0.05说明结果有统计学意义。
plot(anosim_species, col = c('red', "#66C2A5","#F98D63"))

##########标准化###############
###logcmp
library(Seurat)
otu_table_cpm <- NormalizeData(t(phy_species_sub@otu_table),normalization.method = "RC",
                               margin=1,scale.factor = 1000000)
otu_table_cpm <- data.frame(otu_table_cpm,check.names = F)
head(colSums(otu_table_cpm))
otu_table_logcpm <- log(otu_table_cpm+1,2)
head(colSums(otu_table_logcpm))

sample_table <- data.frame(phy_species_sub@sam_data)
sample_table_all <- read.csv("/home/yuekaile/keller/tcga/TCGA-STAD-methylation-phenotype.csv",header = T)

sample_table$gender <- sample_table_all$gender.demographic[match(rownames(sample_table),sample_table_all$submitter_id.samples)]
sample_table$history <- sample_table_all$family_history_of_stomach_cancer[match(rownames(sample_table),sample_table_all$submitter_id.samples)]

sample_table_fil <- sample_table[,c(2,27,33,35,36,29,30)]
colnames(sample_table_fil)[1] <- "age_at_diagnosis"
which(is.na(sample_table_fil),arr.ind = T)
sample_table_fil[c(101,259,313),1] <- round(mean(sample_table_fil$age_at_diagnosis,na.rm=T),0)
sample_table_fil$ostime <- as.numeric(sample_table_fil$ostime)
sample_table_fil$gender <- as.factor(sample_table_fil$gender)
#snm
which(is.na(sample_table_fil),arr.ind = T) 
bio.var <- model.matrix(~os+ostime,
                        data=sample_table_fil)

adj.var <- model.matrix(~age_at_diagnosis +gender+PlateCenter,
                        data=sample_table_fil)

print(dim(adj.var))
print(dim(bio.var))

library(snm)
#####过滤掉和为0的微生物
colsum <- colSums(otu_table_logcpm)
rowsum <- rowSums(otu_table_logcpm)
which(rowsum==0,arr.ind = T)
otu_table_logcpm2 <- otu_table_logcpm %>%
  filter(!rowSums(otu_table_logcpm)==0)
otu_table_logcpm_snm <- snm(as.matrix(otu_table_logcpm2), 
                            bio.var = bio.var, 
                            adj.var = adj.var, 
                            rm.adj=TRUE,
                            verbose = TRUE,
                            diagnose = TRUE)

otu_table_logcpm_snmdata <- data.frame(t(otu_table_logcpm_snm[["norm.dat"]]),check.names = F)
all(rownames(otu_table_logcpm_snmdata)==row.names(sample_table_fil))

data <- cbind(otu_table_logcpm_snmdata,sample_table_fil)

data$history <- ifelse(data$history=="NO","NO",
                       ifelse(data$history=="YES","YES","NO_REPORT"))
which(is.na(data),arr.ind = T) 


###########coxph的单因素回归筛选差异微生物################
cox_result <- data.frame(rep("coef",3424),rep("p",3424))
rownames(cox_result) <- colnames(data)[1:3424]
colnames(cox_result) <- c("coef","p")
data$os2 <- as.numeric(data$os)
for (i in 1:3424) {
  outsur <- coxph(Surv(ostime,os)~data[,i]+data$age_at_diagnosis+data$gender,
                  data=data)
  result <- summary(outsur)
  cox_result[i,1] <- result$coefficients[1,1]
  cox_result[i,2] <- result$coefficients[1,5]
}
cox_result$adj_p <- p.adjust(cox_result$p,method = "fdr")
cox_result_sig <- cox_result[cox_result$adj_p<0.05,]

data2 <- data[,c(1:3424,3430)]
data2$os <- ifelse(data2$os=="0","alive","dead")

data_group <- aggregate(data2[,1:3424],list(data2[,3425]),mean)
data_group <- t(data_group)
colnames(data_group) <- data_group[1,]
data_group <- data.frame(data_group[-1,],check.names = F)

data_group$alive <- as.numeric(data_group$alive)
data_group$dead<- as.numeric(data_group$dead)


data_group$logfc <- log(data_group$dead/data_group$alive,2)
data_group$logfc_abs <- log(abs(data_group$dead/data_group$alive),2)
data_group$delta <- data_group$dead-data_group$alive
cox_result_sig$logfc <- data_group$logfc[match(rownames(cox_result_sig),rownames(data_group))]
cox_result_sig$delta <- data_group$delta[match(rownames(cox_result_sig),rownames(data_group))]
cox_result_sig$logfc_abs <- data_group$logfc_abs[match(rownames(cox_result_sig),rownames(data_group))]
P.Value_t <- 0.5
logFC_t <- 0.5
cox_result_sig$change <- cox_result_sig$logfc_abs

cox_result_sig$change <- ifelse(cox_result_sig$change>0.5,"up",
                                ifelse(cox_result_sig$change<-0.5,"down","not"))
va_plot<- ggplot(data = cox_result_sig,aes(x = logfc_abs,y = -log2(adj_p))) +   
  geom_point(alpha=0.4, size=3.5,aes(color=change)) +   
  ylab("-log2(Pvalue)")+ scale_color_manual(values=c("green", "grey","red"))+   
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +   
  geom_hline(yintercept = -log2(P.Value_t),lty=4,col="black",lwd=0.8) +   
  theme_bw()
va_plot

rownames(cox_result_sig) <- substr(rownames(cox_result_sig),regexpr("s__",rownames(cox_result_sig)),
                                   regexpr("s__",rownames(cox_result_sig))+200)
write.csv(cox_result_sig,file = "cox_result_sig_os_s__.csv")


#############glnment的lasso-cox###############
colnames(data) <- substr(colnames(data),regexpr("s__",colnames(data)),
                                   regexpr("s__",colnames(data))+200)

data_sig <- data[,colnames(data) %in% rownames(cox_result_sig)]

lasso_cox_micro <- as.data.frame(matrix(nrow=100,ncol=101)) 
colnames(lasso_cox_micro) <- rep(1000:1100)
for (i in 1000:1100) {
  set.seed(1046)
train.Idx_sig <- sample(1:dim(data_sig)[1], floor(0.7*dim(data_sig)[1]))
test.Idx_sig <- setdiff(1:dim(data_sig)[1], train.Idx_sig)
x.train_sig <- as.matrix(data_sig[train.Idx_sig ,])
x.test_sig <- as.matrix(data_sig[test.Idx_sig ,])
y.train_sig <- sample_table_fil[train.Idx_sig,]
y.train_y_sig <- data.matrix(Surv(y.train_sig$ostime,y.train_sig$os))
#y.train_x <- y.train
y.test_sig <- sample_table_fil[test.Idx_sig,]
y.test_y_sig <- data.matrix(Surv(y.test_sig$ostime,y.test_sig$os))

fit_sig <- glmnet(x.train_sig,y.train_y_sig,family = "cox",alpha = 1)
plot(fit_sig, xvar = "lambda", label = TRUE)
cvfit_sig = cv.glmnet(x.train_sig,y.train_y_sig, family="cox", alpha=1,nfolds=5) 
plot(cvfit_sig) 
print(cvfit_sig)
coef_sig = coef(fit_sig, s = cvfit_sig$lambda.min) 
index_sig = which(coef_sig != 0) 
actCoef_sig = coef_sig[index_sig] 
lassoGene_sig = row.names(coef_sig)[index_sig] 
if(length(index_sig) > 0){
lasso_cox_micro[1:length(index_sig),i-999] <- lassoGene_sig
}
else{
  next
}
}

# lasso_cox_micro_merge <-gather(lasso_cox_micro)
# lasso_cox_micro_merge <- na.omit(lasso_cox_micro_merge)
# lasso_cox_micro_merge_count <- data.frame(table(lasso_cox_micro_merge$value))
# which(lasso_cox_micro_merge_count$Freq >20)
# lassoGene_sig <- lasso_cox_micro_merge_count[lasso_cox_micro_merge_count$Freq>20,]
# lassoGene_sig_name <- lassoGene_sig$Var1
#actCoef_sig = coef_sig[index_sig]
lassoGene_sig = row.names(coef_sig)[index_sig]
geneCoef_sig = cbind(Gene=lassoGene_sig,Coef=actCoef_sig)
geneCoef_sig
write.csv(data.frame(geneCoef_sig),file = "lasso_cox_stad_os_sig_microlist-1046.csv")


##计算risk score
FinalGeneExp = data_sig[,lassoGene_sig] 
myFun = function(x){crossprod(as.numeric(x),actCoef_sig)}
riskScore = apply(FinalGeneExp,1,myFun)

outCol = c("futime", "fustat", lassoGene_sig) 
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low")) 
data_sig_riskscore = cbind(data_sig[,colnames(data_sig) %in% outCol], riskScore=as.vector(riskScore), risk)

data_sig_riskscore$os <- sample_table_fil$os[match(rownames(data_sig_riskscore),rownames(sample_table_fil))]
library(ggpubr)   #使用ggpubr包绘制散点图
risk_os_plot <- ggboxplot(data_sig_riskscore, x = "os", y = "riskScore", 
                          color = "os", palette = "jco", 
                          add = "jitter") + stat_compare_means()+
  scale_color_manual(values = c("#66C2A5","#F98D63"))
risk_os_plot

#模型的准确性,将lasso筛选出来的cpg作为输入拟合coxph并在验证集上验证
#https://blog.csdn.net/weixin_39954569/article/details/111695267v
#https://blog.csdn.net/weixin_39954569/article/details/111695267
#https://www.jianshu.com/p/60206635270e
library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线 
library(glmnet) 
library(caret)

library(pec) ##验证模型
library(rms)  ##拟合生存分析模型
library(survival)  ##生存分析包
library(glmnet) 
### cox1 为全模型
FinalGeneExp
all(rownames(FinalGeneExp)==rownames(sample_table_fil))
cox_sig13_data <- cbind(FinalGeneExp,sample_table_fil)
cox_sig13_data <- cox_sig13_data[,-c(15,16,18)]

cox_sig13_name <- data.frame(rep(NA,13),rep(NA,13)) 
colnames(cox_sig48_name) <- c("name1","name2")

cox_sig13_name$name1 <- colnames(cox_sig13_data)[1:13]
cox_sig13_name$name2 <- paste("m",rep(1:13),sep = "_")

colnames(cox_sig13_data)[1:13] <- paste("m",rep(1:13),sep = "_")
cox_sig13_data$os <- as.numeric(cox_sig13_data$os)
train_cox <- cox_sig13_data[rownames(x.train_sig),]
test_cox <- cox_sig13_data[rownames(x.test_sig),]
f_cox <- Surv(train_cox$ostime,train_cox$os)

cox_sig_13 <- cph(Surv(ostime,os)~m_1+m_2+m_3+m_4+m_5+m_6+m_7+m_8+m_9+m_10+m_11+m_12+
                  m_13+age_at_diagnosis+gender,data=train_cox, surv=TRUE)

ddist <- datadist(train_cox)
options(datadist='ddist')
surv.cox <- Survival(cox_sig_13)
nom.cox <- nomogram(cox_sig_13, fun=list(function(x) surv.cox(36, x), function(x) surv.cox(60, x), 
                                        function(x) surv.cox(120, x)), lp=F, funlabel=c("3-year survival", "5-year survival", "10-year survival"),
                    maxscale=10, fun.at=c(0.95,0.85,0.75,0.5))
plot(nom.cox)

t <- c(1*365,3*365,5*365)
survprob <- predictSurvProb(cox_sig_13,newd=test_cox,times=t)
head(survprob)
calPolt1 <- calPlot(list("Cox(13 variables)"=cox_sig_13),
                    time=3*365,#设置想要观察的时间点，同理可以绘制其他时间点的曲线
                    data=test_cox,legend.x=0.5,
                    legend.y=0.3,legend.cex=0.8)
print(calPolt1)
calPolt2 <- calPlot(list("Cox(13 variables)"=cox_sig_13),
                    time=5*365,#设置想要观察的时间点，同理可以绘制其他时间点的曲线
                    data=test_cox,legend.x=0.5,
                    legend.y=0.3,legend.cex=0.8)
print(calPolt2)
calPolt3 <- calPlot(list("Cox(13 variables)"=cox_sig_13),
                    time=1*365,#设置想要观察的时间点，同理可以绘制其他时间点的曲线
                    data=test_cox,legend.x=0.5,
                    legend.y=0.3,legend.cex=0.8)
print(calPolt3)


