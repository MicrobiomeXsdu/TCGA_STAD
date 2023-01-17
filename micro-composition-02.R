setwd("/home/yuekaile/keller/stad_rna")
load("~/keller/stad_rna/analysis_dataset.Rdata")

# #############整理成phyloseq######################
# #https://blog.csdn.net/woodcorpse/article/details/106554382
# 
# library(data.table)
# library(phyloseq)
# library(dplyr)
# metadata <- pheno_fil_na
# rownames(metadata) <- metadata$submitter_id_samples
# metadata <- sample_data(metadata)
# otu_genus <- count_g_PlateCenterContamRemoved[-1,]
# otu_genus <- otu_table(as.matrix(otu_genus),taxa_are_rows=TRUE)
# tax_g <- taxnew[taxnew$genus %in% rownames(count_g_PlateCenterContamRemoved),] 
# tax_g <- data.frame(tax_g[,-c(1,8)])
# tax_g <- tax_g[!duplicated(tax_g$genus),]
# rownames(tax_g) <- tax_g$genus
# tax_g <- tax_g[-1,]
# tax_g2 <-  tax_table(as.matrix(tax_g))
# phy_genus <- phyloseq(metadata,otu_genus,tax_g2)
# 
# otu_species <- count_s_PlateCenterContamRemoved[,-1]
# otu_species <- otu_table(as.matrix(otu_species),taxa_are_rows=FALSE)
# 
# tax_s <- taxnew[!taxnew$speices =="s__",] 
# tax_s$name <- substr(tax_s$all,regexpr("g__",tax_s$all),regexpr("g__",tax_s$all)+1000)
# tax_s$name <- gsub("|",".",tax_s$name,fixed = T)
# #tax_s <- data.frame(tax_s[,-1])
# rownames(tax_s) <- tax_s$name
# tax_s2 <- tax_s[rownames(tax_s) %in% colnames(count_s_PlateCenterContamRemoved),]
# 
# tax_s2 <- tax_s2[-1,-8]
# tax_s2 <-  tax_table(as.matrix(tax_s2))
# phy_species <- phyloseq(metadata,otu_species,tax_s2)
# phy_species
# 
# 
# a <- colnames(count_s_PlateCenterContamRemoved)
# a <- data.frame(a)
# b <- a[!duplicated(a$a),]

#############################正式分析--genus水平##########################
#https://www.jianshu.com/p/f7ee51ab69db
load("~/keller/stad_rna/phy_genus_species.Rdata")
###物种组成分析
library(pacman)
library(magrittr)
library(reshape2)
phy_genus@sam_data$groupm <- ifelse(phy_genus@sam_data$pathologic_M=="M1","yes","no")
table(phy_genus@sam_data[["sample_type_samples"]])
# Primary Tumor Solid Tissue Normal 
# 375                  32 
name_match <- read.csv("sample-name.csv",row.names = 1)
name_match$x <- gsub(".","-",name_match$x,fixed = T)
###过滤样本
phy_genus_fill <- prune_samples(phy_genus@sam_data[["submitter_id_samples"]] %in% name_match$x,phy_genus)

table(phy_species_fill@sam_data[["groupm"]])
# no yes 
# 319  19 
table(phy_species_fill@sam_data[["os"]])
# 0   1 
# 195 132 
table(phy_species_fill@sam_data[["groupm"]],phy_species_fill@sam_data[["os"]])
# 0   1
# no  186 123
# yes   9   9
chisq.test(table(phy_species_fill@sam_data[["groupm"]],phy_species_fill@sam_data[["os"]]))
#X-squared = 0.37187, df = 1, p-value = 0.542

pacman::p_load(tidyverse,phyloseq,MicrobiotaProcess,ape)
phytax_genus_l2 <- get_taxadf(obj=phy_genus_fill, taxlevel=2)
phybar_genus_l2 <- ggbartax(obj=phytax_genus_l2,facetNames="groupm", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_genus_l2

phytax_genus_l2m <- get_taxadf(obj=phy_genus_fill, taxlevel=2)
phybar_genus_l2m <- ggbartax(obj=phytax_genus_l2,facetNames="pathologic_M", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_genus_l2m

phytax_genus_l3 <- get_taxadf(obj=phy_genus_fill, taxlevel=3)
phybar_genus_l3 <- ggbartax(obj=phytax_genus_l3,facetNames="groupm", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_genus_l3

phytax_genus_l3m <- get_taxadf(obj=phy_genus_fill, taxlevel=3)
phybar_genus_l3m <- ggbartax(obj=phytax_genus_l3,facetNames="pathologic_M", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_genus_l3m

phytax_genus_l4 <- get_taxadf(obj=phy_genus_fill, taxlevel=4)
phybar_genus_l4 <- ggbartax(obj=phytax_genus_l4,facetNames="groupm", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_genus_l4

phytax_genus_l4m <- get_taxadf(obj=phy_genus_fill, taxlevel=4)
phybar_genus_l4m <- ggbartax(obj=phytax_genus_l4,facetNames="pathologic_M", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_genus_l4m

###########alpha
alpha_genus <- get_alphaindex(phy_genus_fill)
alpha_genus_plot <- ggbox(alpha_genus, geom="violin",factorNames="groupm")+
  scale_fill_manual(values=c("#F7903D","#4D85BD"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
alpha_genus_plot

alpha_genus_plot_m<- ggbox(alpha_genus, geom="violin",factorNames="pathologic_M")+
  scale_fill_manual(values=c("#F7903D","#4D85BD","red"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
alpha_genus_plot_m

alpha_genus_data <- cbind(alpha_genus@alpha,alpha_genus@sampleda)
write.csv(alpha_genus_data,file = "alpha_genus_data.csv")

#######beta多样性+多元置换方差检验
pacman::p_load(tidyverse,ggrepel,vegan,ape,ggsignif,patchwork,multcomp)
otu <-phy_genus_fill@otu_table
groups <- phy_genus_fill@sam_data %>%
  as.list()
#bray距离
t_otu <- t(otu)
pcoa <- vegdist(t(decostand(otu, "norm")),method = "bray") %>% 
  pcoa(correction = "none", rn = NULL)
#unifrac距离
#library(picante)
#pcoa <- picante::unifrac(otu,tree) %>%
#  pcoa(correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
pcoadata <- data.frame(rownames(pcoa$vectors),
                       PC1,PC2,groups$groupm)
colnames(pcoadata) <-c("sample","PC1","PC2","group")
pcoadata$pc1_abs <- abs(pcoadata$PC1)
pcoadata$pc2_abs <- abs(pcoadata$PC2)
pcoadata$pc_sum <- pcoadata$pc1_abs+pcoadata$pc2_abs
write.csv(pcoadata,"pcoadata_groupm-genus.csv")
pcoadata$group <- factor(pcoadata$group,levels=c("no","yes"))
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
  scale_color_manual(values=c("#F7903D","#4D85BD"))+
  stat_ellipse(aes(colour=group),level=0.95,show.legend=F,linetype="dotted")+#aes(colour=group),
  scale_fill_manual(values=c("#F7903D","#4D85BD"))+
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
  scale_fill_manual(values=c("#F7903D","#4D85BD"))+
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
  scale_fill_manual(values=c("#F7903D","#4D85BD"))+
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

otu.adonis=adonis(t_otu~group,data = pcoadata,distance = "bray")
p4 <- ggplot(pcoadata,
             aes(PC1, PC2))+
  geom_text(aes(x = -0.5,
                y = 0.6,
                label = paste("PERMANOVA:\ndf = ",
                              otu.adonis$aov.tab$Df[1],"\nR2 = ",
                              round(otu.adonis$aov.tab$R2[1],4),
                              "\np-value = ",
                              otu.adonis$aov.tab$`Pr(>F)`[1],
                              sep = "")),size = 4) +theme_bw() +
  xlab(NULL) + ylab(NULL) +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
p4
p2+p4+p1+p3 + 
  plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)        

#lesfse找差异菌
library(tidyverse)
library(phyloseq)
library(ggtree)
library(treeio)
library(tidytree)
library(MicrobiotaProcess)
library(psych)
library(pheatmap)
library(reshape2)
set.seed(1001)
diffres_genus_group <- diff_analysis(obj=phy_genus_fill, classgroup="groupm", mlfun="lda",filtermod="pvalue",
                                 firstcomfun = "kruskal.test",firstalpha=0.05,strictmod=F, secondcomfun = "wilcox.test",
                                 subclmin=3,subclwilc=FALSE,secondalpha=0.05, lda=2)
diffres_genus_group

ggdiffbox_genus_group <- ggdiffbox(obj=diffres_genus_group, box_notch=FALSE,factorLevels=list(group=c("no","yes")),colorlist=c("#F7903D","#4D85BD"))
ggdiffbox_genus_group
otu_genus_lefse <- diffres_genus_group@mlres
write.csv(otu_genus_lefse,file = "otu_genus_lefse_name.csv")

#############################正式分析--species水平##########################
#https://www.jianshu.com/p/f7ee51ab69db
load("~/keller/stad_rna/phy_genus_species.Rdata")
###物种组成分析####################
library(pacman)
library(magrittr)
library(reshape2)
phy_species@sam_data$groupm <- ifelse(phy_species@sam_data$pathologic_M=="M1","yes","no")

table(phy_species@sam_data[["sample_type_samples"]])
# Primary Tumor Solid Tissue Normal 
# 375                  32 
name_match <- read.csv("sample-name.csv",row.names = 1)
name_match$x <- gsub(".","-",name_match$x,fixed = T)
phy_species_fill <- prune_samples(phy_species@sam_data[["submitter_id_samples"]] %in% name_match$x,phy_species)



pacman::p_load(tidyverse,phyloseq,MicrobiotaProcess,ape)
phytax_species_l2 <- get_taxadf(obj=phy_species_fill, taxlevel=2)
phybar_species_l2 <- ggbartax(obj=phytax_species_l2,facetNames="groupm", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l2

phytax_species_l2m <- get_taxadf(obj=phy_species_fill, taxlevel=2)
phybar_species_l2m <- ggbartax(obj=phytax_species_l2,facetNames="pathologic_M", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l2m

phytax_species_l3 <- get_taxadf(obj=phy_species_fill, taxlevel=3)
phybar_species_l3 <- ggbartax(obj=phytax_species_l3,facetNames="groupm", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l3

phytax_species_l3m <- get_taxadf(obj=phy_species_fill, taxlevel=3)
phybar_species_l3m <- ggbartax(obj=phytax_species_l3,facetNames="pathologic_M", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l3m

phytax_species_l4 <- get_taxadf(obj=phy_species_fill, taxlevel=4)
phybar_species_l4 <- ggbartax(obj=phytax_species_l4,facetNames="groupm", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l4

phytax_species_l4m <- get_taxadf(obj=phy_species_fill, taxlevel=4)
phybar_species_l4m <- ggbartax(obj=phytax_species_l4,facetNames="pathologic_M", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l4m

phytax_species_l7 <- get_taxadf(obj=phy_species_fill, taxlevel=7)
phybar_species_l7m <- ggbartax(obj=phytax_species_l7,facetNames="groupm", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l7m

phytax_species_l7a <- get_taxadf(obj=phy_species, taxlevel=7)
phybar_species_l7a <- ggbartax(obj=phytax_species_l7a,facetNames="groupm", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")

rela_species_data <- data.frame(phybar_species_l7a[["data"]])
rela_species_data <- rela_species_data[rela_species_data$sample_type_samples=="Primary Tumor",]
rela_species_data2 <- rela_species_data %>%
  group_by(feature) %>%
  summarize(sum = sum(value))

           
           
###########alpha###################################
alpha_species <- get_alphaindex(phy_species_fill)
alpha_species_plot <- ggbox(alpha_species, geom="violin",factorNames="groupm")+
  scale_fill_manual(values=c("#F7903D","#4D85BD"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
alpha_species_plot

alpha_species_plot_m<- ggbox(alpha_species, geom="violin",factorNames="pathologic_M")+
  scale_fill_manual(values=c("#F7903D","#4D85BD","red"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
alpha_species_plot_m

alpha_species_data <- cbind(alpha_species@alpha,alpha_species@sampleda)
write.csv(alpha_species_data,file = "alpha_species_data.csv")

#######beta多样性+adonis\anosim##############################
pacman::p_load(tidyverse,ggrepel,vegan,ape,ggsignif,patchwork,multcomp)
otu <-phy_species_fill@otu_table
groups <- phy_species_fill@sam_data %>%
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
                       PC1,PC2,groups$groupm)
colnames(pcoadata) <-c("sample","PC1","PC2","group")
pcoadata$pc1_abs <- abs(pcoadata$PC1)
pcoadata$pc2_abs <- abs(pcoadata$PC2)
pcoadata$pc_sum <- pcoadata$pc1_abs+pcoadata$pc2_abs
write.csv(pcoadata,"pcoadata_groupm-species.csv")
pcoadata$group <- factor(pcoadata$group,levels=c("no","yes"))
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
  scale_color_manual(values=c("#F7903D","#4D85BD"))+
  stat_ellipse(aes(colour=group),level=0.95,show.legend=F,linetype="dotted")+#aes(colour=group),
  scale_fill_manual(values=c("#F7903D","#4D85BD"))+
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
  scale_fill_manual(values=c("#F7903D","#4D85BD"))+
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
  scale_fill_manual(values=c("#F7903D","#4D85BD"))+
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
                              sep = "")),size = 4) +theme_bw() +
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
anosim_species <-anosim(phy_species_fill@otu_table,phy_species_fill@sam_data[["groupm"]],
                      permutations = 999, distance = "bray")
summary(anosim_species)
# ANOSIM statistic R: 0.0177 
# Significance: 0.238 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0358 0.0483 0.0589 0.0696 
# 
# Dissimilarity ranks between and within classes:
#   0%   25%   50%   75%  100%     N
# Between  7 12588 30133 44740 56913  6061
# no       1 14466 28350 42404 56953 50721
# yes     57  9675 35571 46782 55562   171
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
pairwise.anosim(phy_species_fill@otu_table,phy_species_fill@sam_data[["groupm"]], sim.method="bray", p.adjust.m= "fdr")
# pairs          R p.value p.adjusted
# 1 no vs yes 0.06108336   0.026      0.026
set.seed(555)
pairwise.anosim(phy_species_fill@otu_table,phy_species_fill@sam_data[["pathologic_M"]], sim.method="bray", p.adjust.m= "fdr")
# pairs           R p.value p.adjusted
# 1 M0 vs MX  0.04777519   0.096      0.144
# 2 M0 vs M1  0.06473663   0.027      0.081
# 3 MX vs M1 -0.01857059   0.684      0.684

# * A=0.1881，A>0表示组间差异大于组内差异，A<0表示组间差异小于组内差异。
# * observed delta=0.2667,值越小说明组内差异越小。
# * expected delta=0.3285,值越小说明组间差异越小。
# * significance of delta(P)=0.001，小于0.05说明结果有统计学意义。
plot(anosim_species, col = c('red', "#F7903D","#4D85BD"))


###############lefse找差异菌#########################
library(tidyverse)
library(phyloseq)
library(ggtree)
library(treeio)
library(tidytree)
library(MicrobiotaProcess)
library(psych)
library(pheatmap)
library(reshape2)
set.seed(1001)
diffres_species_group <- diff_analysis(obj=phy_species_fill, classgroup="groupm", mlfun="lda",filtermod="pvalue",
                                     firstcomfun = "kruskal.test",firstalpha=0.05,strictmod=F, secondcomfun = "wilcox.test",
                                     subclmin=3,subclwilc=FALSE,secondalpha=0.05, lda=2)
diffres_species_group

ggdiffbox_species_group <- ggdiffbox(obj=diffres_species_group, box_notch=FALSE,factorLevels=list(group=c("no","yes")),colorlist=c("#F7903D","#4D85BD"))
ggdiffbox_species_group
otu_species_lefse <- diffres_species_group@mlres
write.csv(otu_species_lefse,file = "otu_species_lefse_name.csv")


###########
meta <- data.frame(phy_species_fill@sam_data)
meta <- meta[,-c(1,4,31,32,33)]
meta$os <- as.factor(meta$os)

meta$demo <- ifelse(meta$ethnicity_demographic=="hispanic or latino","hispanic or latino","not hispanic or latino")
library(compareGroups)
res <- compareGroups(groupm ~ ., data = meta)
summary(res)
restab <- createTable(res)
print(restab, which.table = "descr")

all <- descrTable( ~. ,data = meta)
export2xls(all, file='all-stad-rna.xlsx')
export2xls(restab, file='groupm-stad-rna.xlsx')
