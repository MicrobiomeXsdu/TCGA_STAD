setwd("/home/yuekaile/keller/stad_rna")
#############整理stad的metadata######################
library(data.table)
pheno <- fread("TCGA-STAD.GDC_phenotype.tsv")
survi <- fread("TCGA-STAD.survival.tsv")
colnames(pheno)
#提取需要的变量
vari_name <- c("submitter_id.samples","age_at_initial_pathologic_diagnosis",
               "batch_number","submitter_id",
               "days_to_new_tumor_event_after_initial_treatment",
               "lost_follow_up","lymph_node_examined_count",
               "neoplasm_histologic_grade",
               "new_tumor_event_after_initial_treatment",
               "number_of_lymphnodes_positive_by_he",
               "pathologic_M","pathologic_N","pathologic_T",
               "primary_therapy_outcome_success",
               "tissue_source_site","year_of_initial_pathologic_diagnosis",
               "days_to_death.demographic","ethnicity.demographic",
               "ethnicity.demographic","vital_status.demographic",
               "year_of_birth.demographic","site_of_resection_or_biopsy.diagnoses",
               "tumor_stage.diagnoses","disease_type",
               "bcr_id.tissue_source_site","code.tissue_source_site",
               "name.tissue_source_site","sample_type.samples",
               "sample_type_id.samples")
pheno <- data.frame(pheno)
#整合成一个data
pheno_fil <- pheno[,colnames(pheno) %in% vari_name]
pheno_fil$os <- survi$OS[match(pheno_fil$submitter_id.samples,survi$sample)]
pheno_fil$ostime <- survi$OS.time[match(pheno_fil$submitter_id.samples,survi$sample)] 

colnames(pheno_fil) <- gsub(".","_",colnames(pheno_fil),fixed = T)
#读取file-name，并和sample对应
kre_name_stad <- fread("kresult_name-stad.csv")
kre_name_stad <- data.frame(kre_name_stad)
sample_sheet <- fread("gdc_sample_sheet.2022-02-28.tsv")
sample_sheet$file_id <- gsub("_gdc_realn_rehead.bam","",sample_sheet$`File Name`,fixed = T)

sample_sheet<- data.frame(sample_sheet)
colnames(sample_sheet) <- gsub(".","_",colnames(sample_sheet),fixed = T)
kre_name_stad$case_id <- sample_sheet$Sample_ID[match(kre_name_stad$file_id,sample_sheet$file_id)]

pheno_fil$file_name <- kre_name_stad$file_id[match(pheno_fil$submitter_id_samples,kre_name_stad$case_id)]
pheno_fil$file_name[is.na(pheno_fil$file_name)] <- 0
pheno_fil_na <- subset(pheno_fil,!pheno_fil$file_name==0) 

pheno_fil_na$sample_kreken <- kre_name_stad$kresult_id[match(kre_name_stad$case_id,pheno_fil_na$submitter_id_samples)]
write.csv(pheno_fil_na,file = "pheno_fil.csv")

#######################整理otu表############################
otu<- fread("kresult-stad.csv")
rownames(otu) <- otu$V1
tax_raw <- data.frame(otu$V1)
taxnew <- data.frame(otu$V1,"kingdom","phylum","class","order","family","genus","speices")
colnames(taxnew) <- c("all","kingdom","phylum","class","order","family","genus","speices")
for (h in 1:20489) {
  a <- strsplit(tax_raw[h,1],split='|',fixed = T)
  b <- a[[1]]
  cc <- strsplit(b,split=' ',fixed = T)
  for (i in 1:length(cc)) {
    x<- cc[[i]]
    if(substring(x,1,1)=="k"){
      taxnew[h,2]<- x
    }else{
      if(substring(x,1,1)=="p"){
        taxnew[h,3]<- x
      }else{
        if(substring(x,1,1)=="c"){
          taxnew[h,4]<- x
        }else{
          if(substring(x,1,1)=="o"){
            taxnew[h,5]<- x
          }else{
            if(substring(x,1,1)=="f"){
              taxnew[h,6]<- x
            }else{
              if(substring(x,1,1)=="g"){
                taxnew[h,7]<- x
              }else{
                if(substring(x,1,1)=="s"){
                  taxnew[h,8]<- x
              }
              }}}}}}}}
taxnew[taxnew=="kingdom"] <- "k__"
taxnew[taxnew=="phylum"] <- "p__"
taxnew[taxnew=="order"] <- "o__"
taxnew[taxnew=="class"] <- "c__"
taxnew[taxnew=="family"] <- "f__"
taxnew[taxnew=="genus"] <- "g__"
taxnew[taxnew=="speices"] <- "s__"
#write.csv(taxnew,file = "taxonomy_kraken_stad-rna.csv")

otu_tax <- cbind(otu,taxnew)
write.csv(otu_tax,file = "otu_tax_stad-rna.csv")



####################去污的metadata#######################################
library("rjson")
library(dplyr)
# 输入stad的metadata，查看其测序中心,只有一个center：BCGSC
aliquotdata <- fromJSON(file = "metadata.cart.2022-03-01.json")
id <- NA
platform <- NA
center <- NA
file_size <- NA
aliquot_id <- NA
filename <- NA
for (i in 1:407) {
  id[i] <- aliquotdata[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_size[i] <- aliquotdata[[i]][["file_size"]]
  center[i] <- aliquotdata[[i]][["analysis"]][["metadata"]][["read_groups"]][[1]][["sequencing_center"]]
  aliquot_id[i] <- aliquotdata[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  filename[i] <- aliquotdata[[i]][["file_name"]]
}
table(center)
#center:BCGSC :407 

meta_decon <- data.frame(filename, id, file_size,aliquot_id)
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}#从倒数第n为开始取
tmp <- as.character(aliquotdata$aliquot_id)
# 按照plate和测序中心进行污染物筛选
meta_decon$PlateCenter <- factor(substrRight(meta_decon$aliquot_id, 7))

table(meta_decon$PlateCenter)
# 1131-13 1157-13 1602-13 1802-13 1884-13 2055-13 2203-13 2343-13 2402-13 A24K-31 A251-31 A31P-31 
# 10      27      28      25      14      25      12      56      30      13      35      14 
# A32D-31 A33Y-31 A354-31 A36D-31 A39E-31 A414-31 
# 9       7      10      27      14      51 

meta_decon$id_15 <- substr(meta_decon$id,regexpr("TCGA",meta_decon$id),regexpr("TCGA",meta_decon$id)+15)
pheno_fil_na$PlateCenter <- meta_decon$PlateCenter[match(pheno_fil_na$submitter_id_samples,meta_decon$id_15)]
write.csv(pheno_fil_na,file = "pheno_fil.csv")


######################属水平、种水平的otu表###################
count_genus <- fread("kresult_g.csv")
count_genus <- data.frame(count_genus)
rownames(count_genus) <- count_genus$V1
count_genus2 <- count_genus[,-1]
colnames(count_genus2) <- kre_name_stad$case_id

write.csv(count_genus2,file = "count_genus_stad-rna.csv")
count_species <- fread("kresult_s.csv")
count_species2 <- data.frame(count_species)
row.names(count_species2) <- count_species2$V1
count_species2 <- count_species2[,-1]
colnames(count_species2) <- kre_name_stad$case_id
count_species3 <- count_species2 %>%
  summarize(name <- rownames(count_species2),
    sum = rowSums(count_species2)) %>%
  filter(!sum==0)
colnames(count_species3)[1] <- "speices_name"
count_species2 <- count_species2[count_species3$speices_name,]

write.csv(count_species2,file = "count_species_stadrna.csv")


############# 去污 属水平和种水平分开做 ############################## 

# Contamination_simulation_FA.R
# Author: Greg Poore
# Date: Aug 16 2019
# Purpose: To decontaminate TCGA data and see if/how it affects the machine learning analyses

# Load dependencies
require(ggplot2)
require(ggsci)
require(limma)
require(Glimma)
require(edgeR)
require(dplyr)
require(doMC)

numCores <- detectCores()
registerDoMC(cores=numCores)

# 需要的数据：原始的计数数据和aliquot_concentration信息

# NB: SINCE DECONTAM ESSENTIALLY PERFORMS LINEAR REGRESSION BETWEEN READ FRACTIONS AND 
# ANALYTE CONCENTRATIONS, AT LEAST 10 SAMPLES ARE REQUIRED PER PLATE-CENTER COMBINATION
# TO BE PROCESSED FOR IDENTIFYING PUTATIVE CONTAMINANTS. NOTE ALSO THAT ANY CONTAMINANT
# IDENTIFIED IN ANY ONE PLATE-CENTER BATCH WILL BE REMOVED FROM THE WHOLE DATASET
booleanPlateCenter <- as.logical(table(pheno_fil_na$PlateCenter)>10)

sufficientPlateCenter <- names(table(pheno_fil_na$PlateCenter))[booleanPlateCenter]
pheno_fil_na$PlateCenterFlag <- (pheno_fil_na$PlateCenter %in% sufficientPlateCenter)
metadata_obj_PlateCenterSubset <- droplevels(pheno_fil_na[pheno_fil_na$PlateCenterFlag,])

# NB: vbContaminatedDataQC CONTAINS RAW TCGA DATA FROM KRAKEN AND THE SPIKED PSEUDOCONTAMINANTS
name_raw_g <- substr(row.names(count_genus2),regexpr("g__",row.names(count_genus2)),regexpr("g__",row.names(count_genus2))+100)

rownames(count_genus2) <- name_raw_g

count_genus3 <- data.frame(t(count_genus2))
count_genus_PlateCenterSubset <- count_genus3[metadata_obj_PlateCenterSubset$submitter_id_samples,]

# Decontam__
#install.packages("decontam")
require(decontam)
countData <- count_genus_PlateCenterSubset
countMetadata <- metadata_obj_PlateCenterSubset

###############获取由浓度的表格
concer_aliq<- fread("aliquot.tsv")
countMetadata$aliquot_concentration <- concer_aliq$concentration[match(countMetadata$submitter_id_samples,concer_aliq$sample_submitter_id)]

countMetadata$aliquot_concentration <- as.numeric(countMetadata$aliquot_concentration)
countMetadata$aliquot_concentration[is.na(countMetadata$aliquot_concentration)] <- 1000
countMetadata<- countMetadata[!countMetadata$aliquot_concentration==1000,]
rownames(countMetadata) <- countMetadata$submitter_id_samples

countData <- countData[rownames(countMetadata),]
contamdf.freq.plateCenter <- isContaminant(seqtab = as.matrix(countData), 
                                           conc = countMetadata$aliquot_concentration, 
                                           method = "frequency", 
                                           batch = countMetadata$PlateCenter,
                                           threshold = 0.1) # DEFAULT VALUE IS 0.1
save(contamdf.freq.plateCenter, file = "contamdf.freq.plateCenter_stad.RData")
table(contamdf.freq.plateCenter$contaminant)
# FALSE  TRUE 
# 1209    28 


#contamSum <- rownames(contamdf.freq.plateCenter[contamdf.freq.plateCenter$contaminant=="TRUE",])
contamSum <- colSums(as.matrix(count_genus2)[contamdf.freq.plateCenter$contaminant,])
sum(contamSum)/sum(rowSums(as.matrix(count_genus2))) #-->0.000967919
plateCenterSplitGenera <- rownames(contamdf.freq.plateCenter)[contamdf.freq.plateCenter$contaminant]
# generaNamesPlateCenter <- sapply(plateCenterSplitGenera, "[",2)
# tail(generaNamesPlateCenter)
# # NB: This ^ cut off the spiked contaminants, which should be included for downstream processing. FIX WITH FOLLOWING CODE:
# generaNamesPlateCenter[(length(generaNamesPlateCenter)-1):(length(generaNamesPlateCenter))] <-
#   c("contaminant4RandomSpikesHarvard", "contaminant5RandomSpikes1000")

negativa_blank <- c('g__Afipia', 'g__Aquabacteriume', 'g__Asticcacaulis', 'g__Aurantimonas', 'g__Beijerinckia',
                    'g__Bosea', 'g__Bradyrhizobium', 'g__Brevundimonas', 'g__Caulobacter', 'g__Craurococcus', 
                    'g__Devosia', 'g__Hoefleae', 'g__Mesorhizobium', 'g__Methylobacterium', 'g__Novosphingobium', 
                    'g__Ochrobactrum', 'g__Paracoccus', 'g__Pedomicrobium', 'g__Phyllobacteriume',
                    'g__Rhizobium', 'g__Roseomonas', 'g__Sphingobium', 'g__Sphingomonas', 'g__Sphingopyxis',
                    'g__Acidovorax', 'g__Azoarcuse', 'g__Azospira', 'g__Burkholderia', 'g__Comamonas', 'g__Cupriavidus', 
                    'g__Curvibacter', 'g__Delftiae', 'g__Duganella', 'g__Herbaspirillum', 'g__Janthinobacterium', 'g__Kingella', 
                    'g__Leptothrix', 'g__Limnobactere', 'g__Massilia', 'g__Methylophilus', 'g__Methyloversatilise', 
                    'g__Oxalobacter', 'g__Pelomonas', 'g__Polaromonase', 'g__Ralstonia', 'g__Schlegelella', 
                    'g__Sulfuritalea', 'g__Undibacterium', 'g__Variovorax', 'g__Acinetobacter', 'g__Enhydrobacter', 
                    'g__Enterobacter', 'g__Escherichia', 'g__Nevskia', 'g__Pseudomonas', 
                    'g__Pseudoxanthomonas', 'g__Psychrobacter', 'g__Stenotrophomonas', 'g__Xanthomonas', 
                    'g__Aeromicrobium', 'g__Arthrobacter', 'g__Beutenbergia', 'g__Brevibacterium', 'g__Corynebacterium', 
                    'g__Curtobacterium', 'g__Dietzia', 'g__Geodermatophilus', 'g__Janibacter', 'g__Kocuria', 'g__Microbacterium', 
                    'g__Micrococcus', 'g__Microlunatus', 'g__Patulibacter', 'g__Propionibacteriume', 'g__Rhodococcus', 'g__Tsukamurella',
                    'g__Abiotrophia', 'g__Bacillus', 'g__Brevibacillus', 'g__Brochothrix', 'g__Facklamia', 'g__Paenibacillus', 
                    'g__Streptococcus', 'g__Chryseobacterium', 'g__Dyadobacter', 'g__Flavobacterium', 'g__Hydrotalea', 'g__Niastella', 
                    'g__Olivibacter', 'g__Pedobacter', 'g__Wautersiella', 'g__Deinococcus')

popaco <- c('g__Streptococcus', 'g__Mycobacterium', 'g__Staphylococcus', 'g__Waddlia', 'g__Sphingomonas',
            'g__Clostridium', 'g__Porphyromonas', 'g__Micrococcus', 'g__Coxiella', 'g__Fusobacterium', 'g__Yersinia',
            'g__Necropsobacter', 'g__Stenotrophomonas', 'g__Clostridioides', 'g__Borrelia', 'g__Weissella', 'g__Gemella', 'g__Cupriavidus',
            'g__Atopobium', 'g__Brachyspira', 'g__Bifidobacterium', 'g__Lactococcus', 'g__Leptotrichia', 'g__Bartonella', 'g__Circovirus',
            'g__Roseolovirus', 'g__Bilophila', 'g__Borreliella', 'g__Oscillibacter', 'g__Alcaligene', 'g__Ruminococcus',
            'g__Cellulosimicrobium', 'g__Dysgonomonas', 'g__Lachnoanaerobaculum', 'g__Edwardsiella', 'g__Gardnerella',
            'g__Myroides', 'g__Moraxella', 'g__Erysipelothrix', 'g__Plesiomonas', 'g__Odoribacter', 'g__Barnesiella', 'g__Morganella',
            'g__Paraprevotella', 'g__Providencia', 'g__Megasphaera', 'g__Sphingobacterium', 'g__Akkermansia',
            'g__Diplorickettsia', 'g__Empedobacter', 'g__Bergeyella', 'g__Butyricimonas', 'g__Mobiluncus', 'g__Varicellovirus',
            'g__Alloscardovi', 'g__Buttiauxella', 'g__Parachlamydia', 'g__Leclercia', 'g__Gluconobacter',
            'g__Moellerella', 'g__Orientia', 'g__Catabacter', 'g__Cedecea', 'g__Anaerovorax', 'g__Rotavirus', 'g__ Rubivirus')
genus_name <- colnames(count_genus2)
intersect(negativa_blank, popaco)
# "g__Sphingomonas"     "g__Cupriavidus"      "g__Stenotrophomonas"
# "g__Micrococcus"      "g__Streptococcus"
intersect(popaco, rownames(count_genus2)) # 35
intersect(negativa_blank, rownames(count_genus2))# 66
intersect(plateCenterSplitGenera, popaco)# 1
name_inter_blank_popato <- intersect(negativa_blank, popaco)# 5
# RESULTANT PLATE-CENTER DECONTAMINATED DATA:

#####blank-list中不在poti（潜在的致病菌）和stad0.05的genus
blank_list_true <- c("g__Afipia",
                     "g__Aquabacteriume",
                     "g__Bosea",
                     "g__Ochrobactrum",
                     "g__Acidovorax",
                     "g__Curvibacter",
                     "g__Leptothrix",
                     "g__Oxalobacter",
                     "g__Enterobacter",
                     "g__Pseudoxanthomonas",
                     "g__Aeromicrobium",
                     "g__Curtobacterium",
                     "g__Abiotrophia",
                     "g__Olivibacter",
                     "g__Asticcacaulis",
                     "g__Hoefleae",
                     "g__Paracoccus",
                     "g__Azoarcuse",
                     "g__Delftiae",
                     "g__Limnobactere",
                     "g__Pelomonas",
                     "g__Undibacterium",
                     "g__Escherichia",
                     "g__Psychrobacter",
                     "g__Dietzia",
                     "g__Microlunatus",
                     "g__Chryseobacterium",
                     "g__Pedobacter",
                     "g__Aurantimonas",
                     "g__Pedomicrobium",
                     "g__Sphingobium",
                     "g__Azospira",
                     "g__Duganella",
                     "g__Massilia",
                     "g__Polaromonase",
                     "g__Variovorax",
                     "g__Nevskia",
                     "g__Beutenbergia",
                     "g__Geodermatophilus",
                     "g__Patulibacter",
                     "g__Brevibacillus",
                     "g__Dyadobacter",
                     "g__Wautersiella",
                     "g__Beijerinckia",
                     "g__Caulobacter",
                     "g__Phyllobacteriume",
                     "g__Herbaspirillum",
                     "g__Methylophilus",
                     "g__Ralstonia",
                     "g__Xanthomonas",
                     "g__Brevibacterium",
                     "g__Janibacter",
                     "g__Propionibacteriume",
                     "g__Brochothrix",
                     "g__Flavobacterium",
                     "g__Deinococcus",
                     "g__Craurococcus",
                     "g__Sphingopyxis",
                     "g__Comamonas",
                     "g__Janthinobacterium",
                     "g__Methyloversatilise",
                     "g__Schlegelella",
                     "g__Enhydrobacter",
                     "g__Facklamia",
                     "g__Hydrotalea",
                     "g__Kingella",
                     "g__Microbacterium",
                     "g__Tsukamurella",
                     "g__Paenibacillus",
                     "g__Niastella")

#######相当于去掉68个
count_g_PlateCenterContamRemoved <- count_genus2[!(rownames(count_genus2) %in% 
                                                     c(plateCenterSplitGenera, 
                                                       blank_list_true)),]

# # Evaluate how much data was lost
# 1-sum(colSums(count_g_PlateCenterContamRemoved))/sum(colSums(count_genus2))
# # 0.1841489

1-sum(colSums(count_g_PlateCenterContamRemoved[-1,]))/sum(colSums(count_genus2[-1,]))
# 0.1835114

#####只去除找出来的污染物
count_g_PlateCenterContamRemoved1 <- count_genus2[!rownames(count_genus2) %in% 
                                                    plateCenterSplitGenera,]

1-sum(colSums(count_g_PlateCenterContamRemoved1[-1,]))/sum(colSums(count_genus2[-1,]))
# 0.004968305



# 种水平
count_species3 <- count_species2
row.names(count_species3) <- gsub("[|]",".", row.names(count_species3))
name_raw_s <- substr(row.names(count_species3),regexpr("g__",row.names(count_species3)),
                     regexpr("g__",row.names(count_species3))+500)
which(duplicated(name_raw_s))
# integer(0)

# row.names(count_species3)[c(527,530,4522,4523,4524,4525)]
# name_raw_s[c(527,530,4522,4523,4524,4525)]
# # name_raw_s[c(3015, 3016)] <- c("s__Geobacillus_virus_E2", "s__Lactococcus_phage_P335_sensu_lato")

row.names(count_species3) <- name_raw_s

count_species3 <- data.frame(t(count_species3))
count_species_PlateCenterSubset <- count_species3[metadata_obj_PlateCenterSubset$submitter_id_samples,]
# Decontam
require(decontam)
countData <- count_species_PlateCenterSubset
countMetadata <- metadata_obj_PlateCenterSubset
countMetadata$aliquot_concentration <- concer_aliq$concentration[match(countMetadata$submitter_id_samples,concer_aliq$sample_submitter_id)]

countMetadata$aliquot_concentration <- as.numeric(countMetadata$aliquot_concentration)
countMetadata$aliquot_concentration[is.na(countMetadata$aliquot_concentration)] <- 1000
countMetadata<- countMetadata[!countMetadata$aliquot_concentration==1000,]
rownames(countMetadata) <- countMetadata$submitter_id_samples

countData <- countData[rownames(countMetadata),]


contamdf.freq.plateCenter <- isContaminant(seqtab = as.matrix(countData), 
                                           conc = countMetadata$aliquot_concentration, 
                                           method = "frequency", 
                                           batch = countMetadata$PlateCenter,
                                           threshold = 0.1) # DEFAULT VALUE IS 0.1
save(contamdf.freq.plateCenter, file = "contamdf.freq.plateCenter.species.RData")
table(contamdf.freq.plateCenter$contaminant)
# FALSE  TRUE 
# 4457   141 

contamSum <- colSums(as.matrix(count_species3)[,contamdf.freq.plateCenter$contaminant])
sum(contamSum)/sum(colSums(as.matrix(count_species3))) #--> 0.001908856

plateCenterSplitSpecies <- rownames(contamdf.freq.plateCenter)[which(contamdf.freq.plateCenter$contaminant=="TRUE")]
# SpeciesNamesPlateCenter <- sapply(plateCenterSplitSpecies, "[",2)
# tail(SpeciesNamesPlateCenter)
# # NB: This ^ cut off the spiked contaminants, which should be included for downstream processing. FIX WITH FOLLOWING CODE:
# SpeciesNamesPlateCenter[(length(SpeciesNamesPlateCenter)-1):(length(SpeciesNamesPlateCenter))] <- 
#   c("contaminant4RandomSpikesHarvard", "contaminant5RandomSpikes1000")
namespe <- data.frame(t(data.frame(strsplit(name_raw_s, split = ".s__"))))
negativa_blank <- c('g__Afipia', 'g__Aquabacteriume', 'g__Asticcacaulis', 'g__Aurantimonas', 'g__Beijerinckia',
                    'g__Bosea', 'g__Bradyrhizobium', 'g__Brevundimonas', 'g__Caulobacter', 'g__Craurococcus', 
                    'g__Devosia', 'g__Hoefleae', 'g__Mesorhizobium', 'g__Methylobacterium', 'g__Novosphingobium', 
                    'g__Ochrobactrum', 'g__Paracoccus', 'g__Pedomicrobium', 'g__Phyllobacteriume',
                    'g__Rhizobium', 'g__Roseomonas', 'g__Sphingobium', 'g__Sphingomonas', 'g__Sphingopyxis',
                    'g__Acidovorax', 'g__Azoarcuse', 'g__Azospira', 'g__Burkholderia', 'g__Comamonas', 'g__Cupriavidus', 
                    'g__Curvibacter', 'g__Delftiae', 'g__Duganella', 'g__Herbaspirillum', 'g__Janthinobacterium', 'g__Kingella', 
                    'g__Leptothrix', 'g__Limnobactere', 'g__Massilia', 'g__Methylophilus', 'g__Methyloversatilise', 
                    'g__Oxalobacter', 'g__Pelomonas', 'g__Polaromonase', 'g__Ralstonia', 'g__Schlegelella', 
                    'g__Sulfuritalea', 'g__Undibacterium', 'g__Variovorax', 'g__Acinetobacter', 'g__Enhydrobacter', 
                    'g__Enterobacter', 'g__Escherichia', 'g__Nevskia', 'g__Pseudomonas', 
                    'g__Pseudoxanthomonas', 'g__Psychrobacter', 'g__Stenotrophomonas', 'g__Xanthomonas', 
                    'g__Aeromicrobium', 'g__Arthrobacter', 'g__Beutenbergia', 'g__Brevibacterium', 'g__Corynebacterium', 
                    'g__Curtobacterium', 'g__Dietzia', 'g__Geodermatophilus', 'g__Janibacter', 'g__Kocuria', 'g__Microbacterium', 
                    'g__Micrococcus', 'g__Microlunatus', 'g__Patulibacter', 'g__Propionibacteriume', 'g__Rhodococcus', 'g__Tsukamurella',
                    'g__Abiotrophia', 'g__Bacillus', 'g__Brevibacillus', 'g__Brochothrix', 'g__Facklamia', 'g__Paenibacillus', 
                    'g__Streptococcus', 'g__Chryseobacterium', 'g__Dyadobacter', 'g__Flavobacterium', 'g__Hydrotalea', 'g__Niastella', 
                    'g__Olivibacter', 'g__Pedobacter', 'g__Wautersiella', 'g__Deinococcus')

popaco <- c('g__Streptococcus', 'g__Mycobacterium', 'g__Staphylococcus', 'g__Waddlia', 'g__Sphingomonas',
            'g__Clostridium', 'g__Porphyromonas', 'g__Micrococcus', 'g__Coxiella', 'g__Fusobacterium', 'g__Yersinia',
            'g__Necropsobacter', 'g__Stenotrophomonas', 'g__Clostridioides', 'g__Borrelia', 'g__Weissella', 'g__Gemella', 'g__Cupriavidus',
            'g__Atopobium', 'g__Brachyspira', 'g__Bifidobacterium', 'g__Lactococcus', 'g__Leptotrichia', 'g__Bartonella', 'g__Circovirus',
            'g__Roseolovirus', 'g__Bilophila', 'g__Borreliella', 'g__Oscillibacter', 'g__Alcaligene', 'g__Ruminococcus',
            'g__Cellulosimicrobium', 'g__Dysgonomonas', 'g__Lachnoanaerobaculum', 'g__Edwardsiella', 'g__Gardnerella',
            'g__Myroides', 'g__Moraxella', 'g__Erysipelothrix', 'g__Plesiomonas', 'g__Odoribacter', 'g__Barnesiella', 'g__Morganella',
            'g__Paraprevotella', 'g__Providencia', 'g__Megasphaera', 'g__Sphingobacterium', 'g__Akkermansia',
            'g__Diplorickettsia', 'g__Empedobacter', 'g__Bergeyella', 'g__Butyricimonas', 'g__Mobiluncus', 'g__Varicellovirus',
            'g__Alloscardovi', 'g__Buttiauxella', 'g__Parachlamydia', 'g__Leclercia', 'g__Gluconobacter',
            'g__Moellerella', 'g__Orientia', 'g__Catabacter', 'g__Cedecea', 'g__Anaerovorax', 'g__Rotavirus', 'g__ Rubivirus')
negativa_blank1 <- negativa_blank[!negativa_blank%in%popaco]
which(namespe$X1%in%c(negativa_blank1))
which(colnames(count_species)%in%names(contamSum)) #0

intersect(namespe$X1, popaco) # 48
intersect(negativa_blank, rownames(count_species3))

contamSum <- data.frame(contamSum)
# RESULTANT PLATE-CENTER DECONTAMINATED DATA:
#去除593个共计
count_s_PlateCenterContamRemoved <- count_species3[,-unique(c(which(namespe$X1%in%c(blank_list_true)),which(colnames(count_species3)%in%rownames(contamSum))))]

# Evaluate how much data was lost
1-sum(colSums(count_s_PlateCenterContamRemoved[,-1]))/sum(colSums(count_species3[,-1]))
# = 0.1888735

write.csv(count_s_PlateCenterContamRemoved[-1,], "count_s_PlateCenterContamRemoved.csv")
write.csv(count_g_PlateCenterContamRemoved[,-1], "../count_g_PlateCenterContamRemoved.csv")

# saveRDS(c(count_s_PlateCenterContamRemoved[-1,],count_g_PlateCenterContamRemoved[,-1],pheno_fil_na,taxnew,meta_decon),file = "analysis_dataset.RData")
# readRDS("analysis_dataset.RData")

save(count_s_PlateCenterContamRemoved,count_g_PlateCenterContamRemoved,pheno_fil_na,taxnew,meta_decon,file = "analysis_dataset.Rdata")
