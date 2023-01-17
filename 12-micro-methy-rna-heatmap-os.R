dmp_data <- read.csv("/home/yuekaile/keller/ykl/methy-os/coxph-stad-os-methy-dmp-data-12.csv",
                       header = T,check.names = FALSE,row.names = 1)
dmr_data <- read.csv("/home/yuekaile/keller/ykl/methy-os/coxph-stad-os-methy-dmr-data-12.csv",
                     header = T,check.names = FALSE,row.names = 1)
micro_data <- read.csv("/home/yuekaile/keller/stad_rna/micro_os/coxph-stad-os-micro-data-13.csv",
                       header = T,check.names = FALSE,row.names = 1)
micro_data <- micro_data[rownames(dmp_data),]
all(which(rownames(dmp_data)==rownames(dmr_data)))
all(which(rownames(dmp_data)==rownames(micro_data)))
methy_data <- cbind(dmp_data,dmr_data)

library(tidyverse)
library(phyloseq)
library(ggtree)
library(treeio)
library(tidytree)
library(MicrobiotaProcess)
library(psych)
library(pheatmap)
library(reshape2)
###############methy-micro####################################
cor_pre <-corr.test(micro_data[,1:13],dmr_data,method = "spearman",adjust= "none")
cmt_pre <-cor_pre$r
pmt_pre <- cor_pre$p
head(cmt_pre)
head(pmt_pre)
cmt.out_pre<-cbind(rownames(cmt_pre),cmt_pre)
write.table(cmt.out_pre,file= "cor-os-dmr-micro.txt",sep= "t",row.names=F)
pmt.out_pre<-cbind(rownames(pmt_pre),pmt_pre)
write.table(pmt.out_pre,file= "pvalue-os-dmr-micro.txt",sep= "t",row.names=F)
df_pre2 <-reshape2::melt(cmt_pre,value.name= "cor")
df_pre2$pvalue <- as.vector(pmt_pre)
head(df_pre2)
write.table(df_pre2,file= "cor-p-os-dmr-micro.txt",sep= "t")
mycol<-colorRampPalette(c("purple","white","orange"))(200)

if(!is.null(pmt_pre)){
  sssmt <- pmt_pre< 0.01
  pmt_pre[sssmt] <- '**'
  ssmt <- pmt_pre > 0.01& pmt_pre < 0.05
  pmt_pre[ssmt] <- '*'
  smt_pre <- pmt_pre > 0.05& pmt_pre < 0.10
  pmt_pre[smt_pre] <- '.'
  pmt_pre[!sssmt&!ssmt&!smt_pre]<- ''
} else{
  pmt_pre <- F
}
pheat_pre<-pheatmap(cmt_pre,scale = "none",cluster_row = T, cluster_col = T, border=NA,
                    display_numbers = pmt_pre,fontsize_number = 12, number_color = "white",cutoff=0.2,
                    cellwidth = 20, cellheight =20,color=mycol,angle_col = "45",
                     main = "heatmap of corelation between significant microbiota 
                         and immune indice at pre timepoint")


################methy_rna##################
library(ChAMPdata)
probe <- data("probe.features")
library(data.table)
anno <- read.csv("/home/yuekaile/keller/ykl/annotation_gene(1)(1).csv",
                 header = T,check.names = FALSE)
rna_tpm <- fread("/home/yuekaile/keller/ykl/tcga-stad-tpm-allsample.csv",
                 header = T)
rna_tpm <- data.frame(rna_tpm)
row.names(rna_tpm) <- rna_tpm$V1
rna_tpm <- data.frame(t(rna_tpm),check.names = FALSE)
head(rna_tpm[1:5,1:5])
rna_tpm <- rna_tpm[-1,]
rownames(rna_tpm) <- substr(rownames(rna_tpm),1,16)
rownames(rna_tpm) <- gsub(".","-",rownames(rna_tpm),fixed = T)
for (i in 1:56602) {
  rna_tpm[,i] <- as.numeric(rna_tpm[,i])
}
rna_tpm

rna_dmr_gene_name <- unique(probe.features$gene[match(colnames(dmr_data),rownames(probe.features))])
rna_dmr_ensg_name <- anno$gene_id[anno$gene_name %in% rna_dmr_gene_name]
rna_dmr_data <- rna_tpm[rownames(dmr_data),colnames(rna_tpm) %in% rna_dmr_ensg_name]
str(rna_dmr_data)
##################dmr-rna####################
cor_dmr <-corr.test(rna_dmr_data,dmr_data,method = "spearman",adjust= "none")
cmt_dmr <-cor_dmr$r
pmt_dmr <- cor_dmr$p
head(cmt_dmr)
head(pmt_dmr)
cmt.out_dmr<-cbind(rownames(cmt_dmr),cmt_dmr)
write.table(cmt.out_dmr,file= "cor-os-dmr-rna.txt",sep= "t",row.names=F)
pmt.out_dmr<-cbind(rownames(pmt_dmr),pmt_dmr)
write.table(pmt.out_dmr,file= "pvalue-os-dmr-rna.txt",sep= "t",row.names=F)
df_dmr2 <-reshape2::melt(cmt_dmr,value.name= "cor")
df_dmr2$pvalue <- as.vector(pmt_dmr)
head(df_dmr2)
write.table(df_dmr2,file= "cor-p-os-dmr-rna.txt",sep= "t")
mycol<-colorRampPalette(c("purple","white","orange"))(200)

if(!is.null(pmt_dmr)){
  sssmt <- pmt_dmr< 0.01
  pmt_dmr[sssmt] <- '**'
  ssmt <- pmt_dmr > 0.01& pmt_dmr < 0.05
  pmt_dmr[ssmt] <- '*'
  smt_dmr <- pmt_dmr > 0.05& pmt_dmr < 0.10
  pmt_dmr[smt_dmr] <- '.'
  pmt_dmr[!sssmt&!ssmt&!smt_dmr]<- ''
} else{
  pmt_dmr <- F
}
pheat_dmr <-pheatmap(cmt_dmr,scale = "none",cluster_row = T, cluster_col = T, border=NA,
                    display_numbers = pmt_dmr,fontsize_number = 12, number_color = "white",cutoff=0.2,
                    cellwidth = 20, cellheight =20,color=mycol,angle_col = "45",
                    main = "heatmap of corelation between significant microbiota 
                         and immune indice at timepoint")


rna_dmp_gene_name <- unique(probe.features$gene[match(colnames(dmp_data),rownames(probe.features))])
rna_dmp_ensg_name <- anno$gene_id[anno$gene_name %in% rna_dmp_gene_name]
rna_dmp_data <- rna_tpm[rownames(dmp_data),colnames(rna_tpm) %in% rna_dmp_ensg_name]
str(rna_dmp_data)
##################dmp-rna####################
cor_dmp <-corr.test(rna_dmp_data,dmp_data,method = "spearman",adjust= "none")
cmt_dmp <-cor_dmp$r
pmt_dmp <- cor_dmp$p
head(cmt_dmp)
head(pmt_dmp)
cmt.out_dmp<-cbind(rownames(cmt_dmp),cmt_dmp)
write.table(cmt.out_dmp,file= "cor-os-dmp-rna.txt",sep= "t",row.names=F)
pmt.out_dmp<-cbind(rownames(pmt_dmp),pmt_dmp)
write.table(pmt.out_dmp,file= "pvalue-os-dmp-rna.txt",sep= "t",row.names=F)
df_dmp2 <-reshape2::melt(cmt_dmp,value.name= "cor")
df_dmp2$pvalue <- as.vector(pmt_dmp)
head(df_dmp2)
write.table(df_dmp2,file= "cor-p-os-dmp-rna.txt",sep= "t")
mycol<-colorRampPalette(c("purple","white","orange"))(200)

if(!is.null(pmt_dmp)){
  sssmt <- pmt_dmp< 0.01
  pmt_dmp[sssmt] <- '**'
  ssmt <- pmt_dmp > 0.01& pmt_dmp < 0.05
  pmt_dmp[ssmt] <- '*'
  smt_dmp <- pmt_dmp > 0.05& pmt_dmp < 0.10
  pmt_dmp[smt_dmp] <- '.'
  pmt_dmp[!sssmt&!ssmt&!smt_dmp]<- ''
} else{
  pmt_dmp <- F
}
pheat_dmp <-pheatmap(cmt_dmp,scale = "none",cluster_row = T, cluster_col = T, border=NA,
                     display_numbers = pmt_dmp,fontsize_number = 12, number_color = "white",cutoff=0.2,
                     cellwidth = 20, cellheight =20,color=mycol,angle_col = "45",
                     main = "heatmap of corelation between significant microbiota 
                         and immune indice at timepoint")