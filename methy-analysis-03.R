setwd("/home/jinchuandi/ykl")

##################加载环境变量，整理数据###############################
methy_data <- fread("TCGA-STAD.methylation450.tsv")
sample_name <- read.csv("sample-name.csv",header = T,row.names = 1)
sample_name$x <- gsub(".","-",fixed = T,sample_name$x)
colnames(sample_name) <- "name"
methy_data2 <- data.frame(methy_data)
colnames(methy_data2) <- gsub(".","-",colnames(methy_data2),fixed = T)
rownames(methy_data2) <- methy_data2$`Composite-Element-REF`

methy_data2 <- methy_data2[,colnames(methy_data2) %in% sample_name$name]

library(ChAMP)
library("FactoMineR")
library("factoextra")
library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(hgu133plus2.db)
library(limma)
library(org.Hs.eg.db)
library(pheatmap)
require(Biobase)
library("impute")
library(dplyr) 
library(tibble)
#install.packages("pheatmap")
#BiocManager::install("hgu133plus2.db",ask = F,update = F)
#BiocManager::install(c("GSEABase","GSVA"),ask = F,update = F)
#An expression matrix with genes in the rows, samples in the columns
#https://www.jianshu.com/p/1ceace2b581c

#按metadata的顺序排序
methy_data3 <- na.omit(methy_data3)
#485577-375349
#转为矩阵
beta=as.matrix(methy_data3)
#impute missing expression data, k nearest neighbors using a Euclidean metric
beta=impute.knn(beta)
sum(is.na(beta))
#加入伪计数
betaData=beta$data
betaData=betaData+0.00001
a=betaData
a[1:4,1:4]
#表型的metadata

pd <- data.frame(phy_species_fill@sam_data)
pd_2 <- pd[,-c(1,3,4,25,27,28,31)]
#整合成champ的对象
myload_stad = champ.filter(beta = betaData,pd=pd_2)
require(GEOquery)
require(Biobase)
library("impute")
myload_stad[["pd"]][["age_at_initial_pathologic_diagnosis"]] <- as.numeric(myload_stad[["pd"]][["age_at_initial_pathologic_diagnosis"]])
for (i in c(1,2,4,7,13,17,24)) {
  myload_stad[["pd"]][[i]] <- as.numeric(myload_stad[["pd"]][[i]])
}
for (i in c(3,5,6,8:12,14:16,18:23,25:28)) {
  myload_stad[["pd"]][[i]] <- as.factor(myload_stad[["pd"]][[i]])
}

save(myload_stad,file = 'myload_stad_rna.Rdata')


####################标准化前的质控图#########################
library(dplyr) 
library(tibble)
#报错Error in plot.new() : figure margins too large
myload_stad[["pd"]][["groupm"]] <- as.factor(myload_stad[["pd"]][["groupm"]])
QC_pre = champ.QC(beta = myload_stad$beta, pheno = myload_stad$pd$groupm)
#进行标准化
myNorm <- champ.norm(beta=myload_stad$beta,arraytype="450K",cores=10)
#标准化后的pcoa
library(FactoMineR)
library(factoextra)
dat <- t(myNorm)
group_list=myload_stad$pd$groupm
table(group_list)
dat.pca <- PCA(dat, graph = FALSE)
# fviz_pca_ind(dat.pca,geom.ind = "point",col.ind = group_list,addEllipses = TRUE,
#              legend.title = "group")

#检查批次https://www.jianshu.com/p/27269f95d9c1

champ.SVD(beta = myNorm,
          rgSet=NULL,
          pd=myload_stad$pd,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages_batch/")


library(dendextend)#导入R包
dend <- as.dendrogram(hclust(dist(mtcars)))

# #矫正批次效应。
# myCombat <- champ.runCombat(beta=myNorm,pd=myload_stad$pd,variablename=c("groupm","pathologic_N","pathologic_T"),
#                             batchname=c("site_of_resection_or_biopsy_diagnoses"))
# 
# table(myload_stad[["pd"]][["site_of_resection_or_biopsy_diagnoses"]])
# Body of stomach  胃体         Cardia, NOS 胃贲门            Fundus of stomach 胃底
# 79                               84                               43 
# Gastric antrum   胃窦       Lesser curvature of stomach, NOS 胃小弯      Pylorus 幽门
# 122                                1                                        1 
# Stomach, NOS 
# 8 
#write.csv(myCombat,file = "stad_combat_methy_value_data.csv")

##注释到基因
# View(hm450.manifest.hg19)
# View(probe.features)
probe.features <- data.frame(probe.features)
# write.csv(probe.features,file = "450_k_probe.features.csv")

library(dplyr)
library(tidyverse)
probe.features <- separate(data=probe.features,col=UCSC_CpG_Islands_Name,into=c("chorm","start"),sep=":")
probe.features <- separate(data=probe.features,col=start,into=c("start","end"),sep="-")
myCombat_gene <- data.frame(myCombat)
myCombat_gene$gene <- probe.features$gene[match(rownames(myCombat_gene),rownames(probe.features))]


# #查看矫正后 的结果
# champ.SVD(beta = myCombat,
#           rgSet=NULL,
#           pd=myload_stad$pd,
#           RGEffect=FALSE,
#           PDFplot=TRUE,
#           Rplot=TRUE,
#           resultsDir="./CHAMP_SVDimages_after_combat/")



########################差异甲基化位点###########################3
library(tibble)
#b报错myDMP <- champ.DMP(beta = myNorm,pheno=myload_stad$pd$group)
myDMP <- champ.DMP(beta = myNorm,pheno=myload_stad$pd$groupm)
## 4799


# ###加入协变量的差异甲基化位点
# myload_stad[["pd"]][["age_at_initial_pathologic_diagnosis"]] <- as.numeric(myload_stad[["pd"]][["age_at_initial_pathologic_diagnosis"]])
# myload_stad[["pd"]][["gender_class"]] <- ifelse(myload_stad[["pd"]][["gender.demographic"]]=="male",1,2)
# myload_stad[["pd"]][["gender_class"]] <- as.factor(myload_stad[["pd"]][["gender_class"]])
# myDMP_adj <- champ.DMP1(beta =myNorm,pheno = myload_stad$pd$pathologic_M, 
#                         cov1=myload_stad$gender_class,
#                         cov2=myload_stad$age_at_initial_pathologic_diagnosis,
#                         adjust.method = "BH", arraytype = "450K")
# 

head(myDMP[["no_to_yes"]])
df_DMP <- myDMP[["no_to_yes"]] 
df_DMP<- df_DMP[df_DMP$gene!="",]
logFC_t <- 0.05
P.Value_t <- 0.05
summary(df_DMP$deltaBeta)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.13664  0.03972  0.07651  0.07761  0.11291  0.28826 
df_DMP$change <- ifelse(df_DMP$adj.P.Val < P.Value_t & abs(df_DMP$deltaBeta) > logFC_t, 
                        ifelse(df_DMP$deltaBeta > logFC_t ,'UP','DOWN'),'NOT')
table(df_DMP$change) 
# DOWN  NOT   UP 
# 35 1247 2774 
save(df_DMP,file = "df_stad_rna_DMP.Rdata")
#火山图
dataa  = rownames_to_column(df_DMP)
for_label <- dataa%>% head(3)
p <- ggplot(data = dataa,aes(x = logFC,y = -log2(adj.P.Val))) +   
  geom_point(alpha=0.4, size=3.5,aes(color=change)) +   
  ylab("-log2(Pvalue)")+ scale_color_manual(values=c("green", "grey","red"))+   
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +   
  geom_hline(yintercept = -log2(P.Value_t),lty=4,col="black",lwd=0.8) +   
  theme_bw()
p
#
DMP.GUI(DMP=myDMP,beta=myNorm,pheno=myload_stad$pd$groupm)


################差异甲基化区域#################################
myDMR <- champ.DMR(beta=myNorm,pheno=myload_stad$pd$groupm,method="Bumphunter",minProbes=10,cores = 10)
head(myDMR$BumphunterDMR)
df_DMR <- myDMR$BumphunterDMR
write.csv(df_DMR,file = "df_DMR_stad_rna.csv")
df_DMR$names <- as.factor(rownames(df_DMR))
df_DMR$seqnames <- as.factor(df_DMR$seqnames)
DMR.GUI(DMR=myDMR,beta=myNorm,pheno=myload_stad$pd$groupm)
#绘制占比堆砌图

ggplot(df_DMR, aes(x =seqnames , y =L,fill=names)) + 
  geom_col(position = "fill") +
  theme_bw()
# theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
#geom_text(aes(label = L), position = position_stack(vjust = .5), size = 3)  # labels inside the bar segments


##############################差异甲基化block#########################
myBlock <- champ.Block(beta=myNorm,pheno=myload_stad$pd$groupm,arraytype="450K",minNum=100,cores = 10)
head(myBlock$Block)
# library(shinythemes)
# library(plotly)
# Block.GUI(Block=myBlock_batch,beta=myCombat, compare.group=c("M0","M1"), pheno=myload_stad_combat$pd$group,
#           runDMP=TRUE,arraytype="450K")


#########################功能富集分析###################################
myGSEA_fisher <- champ.GSEA(beta=myNorm,DMP=myDMP,DMR=myDMR,CpGlist=NULL,Genelist=NULL,
                     pheno=myload_stad$pd$groupm, method="fisher",arraytype="450K",
                     Rplot=TRUE,adjPval=0.5)
str(myGSEA_fisher)
myGSEA2_goemth <- champ.GSEA(beta=myNorm,DMP=myDMP,DMR=myDMR,CpGlist=NULL,Genelist=NULL,cores=10,
                          pheno=myload_stad$pd$groupm, method="gometh",arraytype="450K",
                          Rplot=TRUE,adjPval=0.5)
str(myGSEA2_goemth)
# barplot(myGSEA_combat_2$DMP, xaxis = "Gene_List", num = 5, colorby = "adjPval", title = "methylglm-KEGG")


###########################差异甲基化模块#############################
##需要跑一下原始脚本224-278.BiocManager::install("FEM")
if(getRversion() >= "2.18.3") utils::globalVariables(c("myNorm","myLoad","hprdAsigH.m"))

champ.EpiMod <- function(beta=myNorm,
                         pheno=myLoad$pd$Sample_Group,
                         nseeds=100,
                         gamma=0.5,
                         nMC=1000,
                         sizeR.v=c(1,100),
                         minsizeOUT=10,
                         resultsDir="./CHAMP_EpiMod/",
                         PDFplot=TRUE,
                         arraytype="450K")
{
  message("[===========================]")
  message("[<<< ChAMP.EpiMod START >>>>]")
  message("-----------------------------")
  
  ### Prepare Checking ###
  if (!file.exists(resultsDir)) dir.create(resultsDir)
  message("champ.EpiMod Results will be saved in ",resultsDir)
  
  data(hprdAsigH)
  message("<< Load PPI network hprdAsigH >>")
  
  if(arraytype=="EPIC")
    statM.o <- GenStatM(beta,pheno,arraytype)
  else
    statM.o <- GenStatM(beta,pheno,"450k")
  message("<< Generate statM.o >>")
  
  intEpi.o=DoIntEpi450k(statM.o,hprdAsigH.m,c=1)
  
  message("<< Calculate EpiMod.o >>")
  EpiMod.o=DoEpiMod(intEpi.o,
                    nseeds=nseeds,
                    gamma=gamma,
                    nMC=nMC,
                    sizeR.v=sizeR.v,
                    minsizeOUT=minsizeOUT,
                    writeOUT=TRUE,
                    ew.v=NULL);
  
  if(PDFplot)
  {
    message("<< Draw All top significant module plot in PDF >>")
    tmpdir <- getwd()
    setwd(resultsDir)
    for(i in names(EpiMod.o$topmod)) FemModShow(EpiMod.o$topmod[[i]],name=i,EpiMod.o,mode="Epi")
    setwd(tmpdir)
  }
  
  message("[<<<< ChAMP.EpiMod END >>>>>]")
  message("[===========================]")
  return(EpiMod.o)
}

library(FEM)
data(hprdAsigH)
class(hprdAsigH.m)
myepimod <- champ.EpiMod(beta=myNorm,pheno=myload_stad$pd$groupm,nseeds=100,
                         gamma=0.5,nMC=1000,sizeR.v=c(1,100),minsizeOUT=10,
                         resultsDir="./CHAMP_EpiMod/",PDFplot=TRUE,arraytype="450K")
df_epimod <- myepimod$fem


#####################整体的甲基化水平####################################
total <- data.frame(t(myload_stad[["beta"]]))
total$all <- rowSums(total)/(ncol(total))
total$group <- pd$groupm[match(rownames(total),rownames(pd))]
tapply(total$all,total$group,shapiro.test)
library(ggpubr)
all_methy_level <- ggplot(total, aes(group,all)) +
  geom_point(aes(colour=group),size=2.5)+
  geom_boxplot(aes(colour=group),size=2.5)+
  scale_color_manual(values=c("#228B22","#FF4500"))+
  stat_compare_means()
all_methy_level
write.csv(total_combat2,file = "stad_total_methydata.csv")

#
#save(myload_stad_combat,myDMP3,myDMR2,myGSEA_combat_2,myepimod_batch,myBlock_batch,myCombat,file = "ChAMP_DMP_DMR_GSEA_EPIMOD_version2_combat_0.10.Rdata")
#将mod中的基因表整体出来
mod_gene <- read.table("topEPI-Epi-X.txt",sep = "\t",header = T)
total_methy <- data.frame(myload_stad[["beta"]])
total_methy$gene <- probe.features$gene[match(rownames(total_methy),rownames(probe.features))]

name_csf1 <-c("CSF1","CLSTN1","TXN","TNR","IL16","NPTX2",
              "CLEC11A","TMPRSS11A","LCN12","ADAMTS20",
              "XYLT1","CSF1R","SLK","APBA2","NPTXR")
mod_csf1 <- total_methy[total_methy$gene %in% name_csf1,]

name_nfe2 <- c("NFE2","PEX14","NFE2L2","MAFK","PEX7","BACH2","PATZ1",
               "BACH1","PEX5","PEX12","NFE2L3","PEX10","PEX19","PXMP4",
               "SLC25A17","ABCD3","PEX16","PEX26","PEX6","PEX2")
mod_nfe2  <- total_methy[total_methy$gene %in% name_nfe2,]

name_cenpf <- c("CENPF","PHC2","FER","KHDRBS2","PHC1","ERCC6","MDM1","ZNF546",
                  "NUP133","TOP3B","FNTB","SFMBT1","ERCC5","XPA","TMEM70")
mod_cenpf  <- total_methy[total_methy$gene %in% name_cenpf,]

name_NOTCH1  <- c("NOTCH1","CCN3","HEY2","NUMBL","MAML1","PSEN1","PSEN2",
"ADAM17","APBA1","MFAP5","JAG1","DLL1","APH1A","CNTN1",
"APH1B","JAG2","DNER","CNTN6","MAML3","MFNG","ADAM10","LFNG",
"DTX1","DLK1","PTCRA","NOTCH3","NOTCH4","APP","ICAM5","KCNIP3",
"DOCK3","MEGF6","CNTNAP4","NEURL1","KCNAB1","CD177","MTCH1","GATA6",
"DTX3","KIF17","LIN7A","TMED10","CTNND2","MAGI2","SLC39A2","TMED9",
"DSCAML1","PLCL2")
mod_NOTCH1  <- total_methy[total_methy$gene %in% name_NOTCH1,]

name_RAP2A <- c("RAP2A","RAPGEF6","RASSF5","TNIK","ARHGAP29","RUNDC3A","MAP4K4",
"RGS14","FNTB","RAPGEF5","MRAS","KIF26B","RAP2B")
mod_RAP2A  <- total_methy[total_methy$gene %in% name_RAP2A,]

name_PSEN2 <- c("PSEN2","NOTCH3","PSEN1","APP","NOTCH4","CTNND2","APH1A","APH1B",
          "JAG2","DOCK3","KCNIP3","ICAM5","DNER","CNTN6","DLL1","MFNG","CD177",
          "MTCH1","MAGI2","CNTN1","KCNAB1","TMED10","SLC39A2","TMED9","DSCAML1",
          "PLCL2")
mod_PSEN2  <- total_methy[total_methy$gene %in% name_PSEN2,]

name_FEZ1 <- c("FEZ1","BRD1","PTN","P4HB","TAC3","NBR1","COL9A2","WWC1","PDCD7",
"NEK1","PTH","CLASP2","RAB3GAP2","OLFML3","PAM16","STAR","SCOC","PPP4R3B",
"RAB3GAP1","DDN","WFDC2","CLIP2","KIAA0513","CTNNAL1","TACR2","TAC4","TACR3",
"TACR1","TAC1","INTS4","ZNF350","PTH2R","PTH2")
mod_FEZ1  <- total_methy[total_methy$gene %in% name_FEZ1,]

write.csv(mod_csf1,file = "stad_rna_mod_csf1.csv")
write.csv(mod_nfe2,file = "stad_rna_mod_nfe2.csv")
write.csv(mod_cenpf,file = "stad_rna_mod_cenpf.csv")
write.csv(mod_NOTCH1,file = "stad_rna_mod_NOTCH1.csv")
write.csv(mod_RAP2A,file = "stad_rna_mod_RAP2A.csv")
write.csv(mod_PSEN2,file = "stad_rna_mod_PSEN2.csv")
write.csv(mod_FEZ1,file = "stad_rna_mod_FEZ1.csv")

save(methy_data3,myload_stad,pd,myNorm,probe.features,
     myDMP,df_DMP,myDMR,df_DMR,myBlock,myGSEA_fisher,
     myGSEA2_goemth,myepimod,df_epimod,total,mod_gene,
     mod_csf1,mod_nfe2,mod_cenpf,mod_NOTCH1,mod_RAP2A,
     mod_PSEN2,mod_FEZ1,file = "stad_rna_main_output.Rdata")
