setwd("/home/yuekaile/keller/stad_rna/tumor-normal-result/")
table(phy_species@sam_data[["sample_type_samples"]])
load("~/keller/stad_rna/phy_genus_species.Rdata")
# Primary Tumor Solid Tissue Normal 
# 375                  32 


pacman::p_load(tidyverse,phyloseq,MicrobiotaProcess,ape)
phytax_genus_l2 <- get_taxadf(obj=phy_genus, taxlevel=2)
phybar_genus_l2 <- ggbartax(obj=phytax_genus_l2,facetNames="sample_type_samples", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_genus_l2

phytax_genus_l3 <- get_taxadf(obj=phy_genus, taxlevel=3)
phybar_genus_l3 <- ggbartax(obj=phytax_genus_l3,facetNames="sample_type_samples", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_genus_l3

phytax_genus_l4 <- get_taxadf(obj=phy_genus, taxlevel=4)
phybar_genus_l4 <- ggbartax(obj=phytax_genus_l4,facetNames="sample_type_samples", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_genus_l4

###########alpha
alpha_genus <- get_alphaindex(phy_genus)
alpha_genus_plot <- ggbox(alpha_genus, geom="violin",factorNames="sample_type_samples")+
  scale_fill_manual(values=c("#909FCA","#ACD851"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
alpha_genus_plot

alpha_genus_data <- cbind(alpha_genus@alpha,alpha_genus@sampleda)
write.csv(alpha_genus_data,file = "alpha_genus_data_tumor-normal.csv")

#######beta多样性+多元置换方差检验
pacman::p_load(tidyverse,ggrepel,vegan,ape,ggsignif,patchwork,multcomp)
otu <-phy_genus@otu_table
groups <- phy_genus@sam_data %>%
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
                       PC1,PC2,groups$sample_type_samples)
colnames(pcoadata) <-c("sample","PC1","PC2","group")
pcoadata$pc1_abs <- abs(pcoadata$PC1)
pcoadata$pc2_abs <- abs(pcoadata$PC2)
pcoadata$pc_sum <- pcoadata$pc1_abs+pcoadata$pc2_abs
write.csv(pcoadata,"pcoadata_genus_tvsm.csv")
pcoadata$group <- factor(pcoadata$group)
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
  scale_color_manual(values=c("#909FCA","#ACD851"))+
  stat_ellipse(aes(colour=group),level=0.95,show.legend=F,linetype="dotted")+#aes(colour=group),
  scale_fill_manual(values=c("#909FCA","#ACD851"))+
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
  scale_fill_manual(values=c("#909FCA","#ACD851"))+
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
  scale_fill_manual(values=c("#909FCA","#ACD851"))+
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
                              "\np = ",
                              otu.adonis$aov.tab$`Pr(>F)`[1],
                              sep = "")),size =2.5) +theme_bw() +
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
diffres_genus_group <- diff_analysis(obj=phy_genus, classgroup="sample_type_samples", mlfun="lda",filtermod="pvalue",
                                     firstcomfun = "kruskal.test",firstalpha=0.05,strictmod=F, secondcomfun = "wilcox.test",
                                     subclmin=3,subclwilc=FALSE,secondalpha=0.05, lda=3)
diffres_genus_group

ggdiffbox_genus_group <- ggdiffbox(obj=diffres_genus_group, box_notch=FALSE,
                                   colorlist=c("#909FCA","#ACD851"))
ggdiffbox_genus_group
otu_genus_lefse <- diffres_genus_group@mlres
write.csv(otu_genus_lefse,file = "otu_genus_lefse_name_tvsn.csv")

#############################正式分析--species水平##########################
#https://www.jianshu.com/p/f7ee51ab69db
load("~/keller/stad_rna/phy_genus_species.Rdata")
###物种组成分析####################
library(pacman)
library(magrittr)
library(reshape2)
table(phy_species@sam_data[["sample_type_samples"]])
# Primary Tumor Solid Tissue Normal 
# 375                  32 

pacman::p_load(tidyverse,phyloseq,MicrobiotaProcess,ape)
phytax_species_l2 <- get_taxadf(obj=phy_species, taxlevel=2)
phybar_species_l2 <- ggbartax(obj=phytax_species_l2,facetNames="sample_type_samples", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l2

phytax_species_l3 <- get_taxadf(obj=phy_species, taxlevel=3)
phybar_species_l3 <- ggbartax(obj=phytax_species_l3,facetNames="sample_type_samples", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l3


phytax_species_l4 <- get_taxadf(obj=phy_species, taxlevel=4)
phybar_species_l4 <- ggbartax(obj=phytax_species_l4,facetNames="sample_type_samples", count=FALSE) +
  xlab(NULL) + ylab("relative abundance (%)")+
  theme(axis.text.x=element_text(face="plain",
                                 color="black",hjust=0.8,vjust=0.6,
                                 size=9, angle=90))+
  theme(strip.text.x = element_text(size=8, color="black",
                                    face="plain"))+
  theme(legend.position="right")
phybar_species_l4

###########alpha###################################
alpha_species <- get_alphaindex(phy_species)
alpha_species_plot <- ggbox(alpha_species, geom="violin",factorNames="sample_type_samples")+
  scale_fill_manual(values=c("#909FCA","#ACD851"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
alpha_species_plot
alpha_species_data <- cbind(alpha_species@alpha,alpha_species@sampleda)
write.csv(alpha_species_data,file = "alpha_species_data_tvsn.csv")

#######beta多样性+adonis\anosim##############################
pacman::p_load(tidyverse,ggrepel,vegan,ape,ggsignif,patchwork,multcomp)
otu <-phy_species@otu_table
groups <- phy_species@sam_data %>%
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
                       PC1,PC2,groups$sample_type_samples)
colnames(pcoadata) <-c("sample","PC1","PC2","group")
pcoadata$pc1_abs <- abs(pcoadata$PC1)
pcoadata$pc2_abs <- abs(pcoadata$PC2)
pcoadata$pc_sum <- pcoadata$pc1_abs+pcoadata$pc2_abs
write.csv(pcoadata,"pcoadata_sample_type_samples-species.csv")
pcoadata$group <- factor(pcoadata$group)
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
  scale_color_manual(values=c("#909FCA","#ACD851"))+
  stat_ellipse(aes(colour=group),level=0.95,show.legend=F,linetype="dotted")+#aes(colour=group),
  scale_fill_manual(values=c("#909FCA","#ACD851"))+
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
  scale_fill_manual(values=c("#909FCA","#ACD851"))+
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
  scale_fill_manual(values=c("#909FCA","#ACD851"))+
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
                              "\np = ",
                              otu.adonis$aov.tab$`Pr(>F)`[1],
                              sep = "")),size = 2.5) +theme_bw() +
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
anosim_species <-anosim(phy_species@otu_table,phy_species@sam_data[["sample_type_samples"]],
                        permutations = 999, distance = "bray")
summary(anosim_species)
# Call:
#   anosim(x = phy_species@otu_table, grouping = phy_species@sam_data[["sample_type_samples"]],      permutations = 999, distance = "bray") 
# Dissimilarity: bray 
# 
# ANOSIM statistic R: -0.1042 
# Significance: 1 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0402 0.0506 0.0619 0.0743 
# 
# Dissimilarity ranks between and within classes:
#   0%     25%     50%      75%  100%     N
# Between              13 19562.5 35538.5 55915.00 82405 12000
# Primary Tumor         1 20980.0 42600.0 63118.00 82621 70125
# Solid Tissue Normal 514 15365.5 26048.0 40434.25 72199   496

# dc_species = data.frame(dis=anosim_species$dis.rank, class=anosim_species$class.vec)
# # 多重比较
# pairwise.anosim <-function(x,factors, sim.method, p.adjust.m){
#   library(vegan)
#   co = as.matrix(combn(unique(factors),2))
#   pairs = c()
#   R = c()
#   p.value = c()  
#   for(elem in 1:ncol(co)){
#     ad = anosim(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),],                
#                 factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))], permutations = 999,distance = "bray");    
#     pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));    
#     R = c(R,ad$statistic);   
#     p.value = c(p.value,ad$signif)  
#   }
#   p.adjusted =p.adjust(p.value,method=p.adjust.m) 
#   pairw.res = data.frame(pairs, R, p.value,p.adjusted) 
#   return(pairw.res)
# } # 不需要更改
# set.seed(5555)
# pairwise.anosim(phy_species@otu_table,phy_species@sam_data[["sample_type_samples"]], sim.method="bray", p.adjust.m= "fdr")
# # pairs          R p.value p.adjusted
# # 1 no vs yes 0.06108336   0.026      0.026
# set.seed(555)
# pairwise.anosim(phy_species@otu_table,phy_species@sam_data[["pathologic_M"]], sim.method="bray", p.adjust.m= "fdr")
# # pairs           R p.value p.adjusted
# # 1 M0 vs MX  0.04777519   0.096      0.144
# # 2 M0 vs M1  0.06473663   0.027      0.081
# # 3 MX vs M1 -0.01857059   0.684      0.684
# 
# # * A=0.1881，A>0表示组间差异大于组内差异，A<0表示组间差异小于组内差异。
# # * observed delta=0.2667,值越小说明组内差异越小。
# # * expected delta=0.3285,值越小说明组间差异越小。
# # * significance of delta(P)=0.001，小于0.05说明结果有统计学意义。
plot(anosim_species, col = c("#909FCA","#ACD851"))


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
diffres_species_group <- diff_analysis(obj=phy_species, classgroup="sample_type_samples", mlfun="lda",filtermod="pvalue",
                                       firstcomfun = "kruskal.test",firstalpha=0.05,strictmod=F, secondcomfun = "wilcox.test",
                                       subclmin=3,subclwilc=FALSE,secondalpha=0.05, lda=3)
diffres_species_group

otu_species_lefse <- diffres_species_group@mlres$f
otu_species_lefse_g <-  otu_species_lefse[17:27]

ggdiffbox_species_group_g <- ggdiffbox(obj=diffres_species_group, box_notch=FALSE,
                                       featurelist=otu_species_lefse_g,
                                     colorlist=c("#909FCA","#ACD851"))
ggdiffbox_species_group_g

otu_species_lefse_s <-  otu_species_lefse[36:49]

ggdiffbox_species_group_s <- ggdiffbox(obj=diffres_species_group, box_notch=FALSE,
                                       featurelist=otu_species_lefse_s,
                                       colorlist=c("#909FCA","#ACD851"))
ggdiffbox_species_group_s

write.csv(otu_species_lefse,file = "otu_species_lefse_name_tvsn.csv")


diffres_species_group_2 <- diff_analysis(obj=phy_species, classgroup="sample_type_samples", mlfun="lda",filtermod="pvalue",
                                       firstcomfun = "kruskal.test",firstalpha=0.05,strictmod=F, secondcomfun = "wilcox.test",
                                       subclmin=3,subclwilc=FALSE,secondalpha=0.05, lda=2)
diffres_species_group_2

otu_species_lefse2 <- diffres_species_group_2@mlres$f
otu_species_lefse2

otu_species_lefse_g22 <-  otu_species_lefse2[c(44:81,177)]

ggdiffbox_species_group_g_2 <- ggdiffbox(obj=diffres_species_group_2, box_notch=FALSE,
                                       featurelist=otu_species_lefse_g22,
                                       colorlist=c("#909FCA","#ACD851"))
ggdiffbox_species_group_g_2

otu_species_lefse_s22 <-  otu_species_lefse2[c(106:176,178:186)]

ggdiffbox_species_group_s2 <- ggdiffbox(obj=diffres_species_group_2, box_notch=FALSE,
                                       featurelist=otu_species_lefse_s22,
                                       colorlist=c("#909FCA","#ACD851"))
ggdiffbox_species_group_s2


otu_data_407_normal <- data.frame(diffres_species_group_2@originalD)
write.csv(otu_data_407_normal,file="otu_data_407_normal.csv")

write.csv(otu_species_lefse2,file = "otu_species_lefse_name_tvsn_05-2.csv")

###########summary description######
meta <- data.frame(phy_species@sam_data)
meta <- meta[,-c(1,4,31,32,33)]
meta$os <- as.factor(meta$os)

meta$demo <- ifelse(meta$ethnicity_demographic=="hispanic or latino","hispanic or latino","not hispanic or latino")
library(compareGroups)
res <- compareGroups(sample_type_samples ~ ., data = meta)
summary(res)
restab <- createTable(res)
print(restab, which.table = "descr")

all <- descrTable( ~. ,data = meta)
export2xls(all, file='all-stad-rna_tumor_normal.xlsx')
export2xls(restab, file='sample_type_samples-stad-rna_tumor_normal.xlsx')


