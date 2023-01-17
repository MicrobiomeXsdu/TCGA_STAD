##https://blog.csdn.net/dingming001/article/details/72909836?spm=1001.2101.3001.6650.14&utm_medium=distribute.pc_relevant.none-task-blog-2~default~BlogCommendFromBaidu~Rate-14.pc_relevant_default&depth_1-utm_source=distribute.pc_relevant.none-task-blog-2~default~BlogCommendFromBaidu~Rate-14.pc_relevant_default&utm_relevant_index=20
setwd("/home/yuekaile/keller/stad_rna/model/")
load("~/keller/stad_rna/phy_genus_species.Rdata")
name <- read.csv("/home/yuekaile/keller/stad_rna/tumor-normal-result/otu_species_lefse_name_tvsn_05-2.csv",header = T,row.names = 1)
library(stringr)
species_name <- name[grep(pattern="s__",name[,1]),]

otu_table <- data.frame(phy_species@otu_table)
colnames(otu_table) <- substr(colnames(otu_table),regexpr("s__",colnames(otu_table)),
                              regexpr("s__",colnames(otu_table))+300)

###logcmp
library(Seurat)
otu_table_cpm <- NormalizeData(t(otu_table),normalization.method = "RC",
                                  margin=1,scale.factor = 1000000)
otu_table_cpm <- data.frame(otu_table_cpm,check.names = F)
head(colSums(otu_table_cpm))
otu_table_logcpm <- log(otu_table_cpm+1,2)
head(colSums(otu_table_logcpm))

sample_table <- data.frame(phy_species@sam_data)

sample_table_all <- read.csv("/home/yuekaile/keller/tcga/TCGA-STAD-methylation-phenotype.csv",header = T)

sample_table$gender <- sample_table_all$gender.demographic[match(rownames(sample_table),sample_table_all$submitter_id.samples)]
sample_table$history <- sample_table_all$family_history_of_stomach_cancer[match(rownames(sample_table),sample_table_all$submitter_id.samples)]

sample_table_fil <- sample_table[,c(2,27,33,35,36)]
colnames(sample_table_fil)[1] <- "age_at_diagnosis"
sample_table_fil[c(120,303,362,377),1] <- 66

#snm
which(is.na(sample_table_fil),arr.ind = T) 
bio.var <- model.matrix(~sample_type_samples,
                        data=sample_table_fil)

adj.var <- model.matrix(~age_at_diagnosis +gender+PlateCenter,
                        data=sample_table_fil)

print(dim(adj.var))
print(dim(bio.var))

library(snm)
otu_table_logcpm_snm <- snm(as.matrix(otu_table_logcpm), 
                      bio.var = bio.var, 
                      adj.var = adj.var, 
                      rm.adj=TRUE,
                      verbose = TRUE,
                      diagnose = TRUE)

otu_table_logcpm_snmdata <- data.frame(t(otu_table_logcpm_snm[["norm.dat"]]),check.names = F)

otu_sig <- otu_table_logcpm_snmdata[,colnames(otu_table_logcpm_snmdata) %in% species_name]
all(rownames(otu_sig)==row.names(sample_table_fil))

data <- cbind(otu_sig,sample_table_fil)

data$history <- ifelse(data$history=="NO","NO",
                       ifelse(data$history=="YES","YES","NO_REPORT"))
which(is.na(data),arr.ind = T) 

###################建模##########
library(lattice)
library(ggplot2)
library(caret)
#churn <- read.csv("/home/yuekaile/keller/stad_rna/model/churn.csv",header = T)
# str(churnTrain)
# churnTrain = churnTrain[,!names(churnTrain) %in% c("state","area_code","account_length")]
# #生成随机编号为2的随机数
# set.seed(2222)
# #将churnTrain的数据集分为两类，按0.7与0.3的比例无放回抽样
# ind = sample(2,nrow(data),replace = TRUE,prob = c(0.7,0.3))
# 
# trainset = data[ind == 1,]
# testset = data[ind == 2,]
# 
# control = trainControl(method = "repeatedcv",number = 10,repeats = 3)
# library(rpart)
# library(C50)
# 
# model = train(sample_type_samples~.,data = trainset,method = "rpart",
#               preProcess = "scale" ,trControl = control)
# 
# importance = varImp(model,scale = FALSE)
# importance


###############特征选择————筛选important的菌###########
#4_feature_selection
#BiocManager::install("optparse")
library(optparse)
library(dplyr)
library(ggplot2)
library(randomForest)
library(caret)
library(A3)
library(Boruta)
library(ROCR)
#library(ImageGP)
library(corrplot)
library(vegan)
library(Hmisc)
###CHECK PACKAGES
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# for(p in package_list){
#   if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
#     install.packages(p, repos=site)
#     suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
#   }
# }

###feature_selection
set.seed(1)
data$sample_type_samples <- gsub(" ","_",data$sample_type_samples,fixed = T)
data$sample_type_samples <- as.factor(data$sample_type_samples)

boruta <- Boruta(sample_type_samples~.,data = data , pValue=0.05, mcAdj=T,
                 maxRuns=1000)
boruta
print(table(boruta$finalDecision))
#extract feature
plot(boruta)
boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = F), 
                                            Type="Boruta_with_tentative")
feature <- boruta.finalVarsWithTentative$Item
boruta.variable.imp.use <- boruta$ImpHistory[,feature]
feature_importance <- apply(boruta.variable.imp.use,2,mean)
feature_importance <- data.frame(sort(feature_importance,decreasing = TRUE))
feature <- rownames(feature_importance)
write.csv(feature_importance,'lefse_species_s__boruta_feature_importance_normalization.csv')

##importance
boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision),
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>%
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}
boruta.variable.imp <- boruta.imp(boruta)
TentativeRoughFix(boruta)
# Boruta performed 308 iterations in 15.25427 secs.
# 16 attributes confirmed important: s__Acinetobacter_sp._WCHA45,
# s__Aliarcobacter_cryaerophilus, s__Brachybacterium_faecium,
# s__Brachybacterium_ginsengisoli, s__Brachybacterium_saurashtrense
# and 11 more;
# 68 attributes confirmed unimportant: age_at_diagnosis, gender,
# history, PlateCenter, s__Acinetobacter_baumannii and 63 more;

getConfirmedFormula(boruta)
# sample_type_samples ~ s__Acinetobacter_sp._WCHA45 + s__Sphingomonas_melonis + 
#   s__Celeribacter_indicus + s__Haematobacter_massiliensis + 
#   s__Aliarcobacter_cryaerophilus + s__Helicobacter_pylori + 
#   s__Staphylococcus_cohnii + s__Brachybacterium_sp._SGAir0954 + 
#   s__Brachybacterium_saurashtrense + s__Brachybacterium_faecium + 
#   s__Brachybacterium_ginsengisoli + s__Brachybacterium_vulturis + 
#   s__Micrococcus_luteus + s__Corynebacterium_ureicelerivorans + 
#   s__Corynebacterium_ammoniagenes + s__Human_mastadenovirus_C  
  
data_boruta <- data[,colnames(data) %in% boruta.finalVarsWithTentative$Item]
data_boruta <- cbind(data_boruta,data$sample_type_samples)
colnames(data_boruta)[17] <- "Group"
#data_boruta$Group <- gsub(" ","_",data_boruta$Group,fixed = T)

##########classification_model_cv_repeat(05)####
# define customRF
customRF <- list(type = "Classification",
                 library = "randomForest",
                 loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree","nodesize","maxnodes"),
                                  class = rep("numeric", 4),
                                  label = c("mtry", "ntree","nodesize","maxnodes"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

### cv fold
set.seed(123)
tune_model <- function(data,Group,num){
  print("tune model")
  set.seed(2021*5)
  control <- trainControl(method="repeatedcv", number = 5,classProbs=T,summaryFunction=twoClassSummary)
  tunegrid <- expand.grid(.mtry=c(1:ceiling(sqrt(length(feature)))),
                          .ntree=c(301,501,801,1001),
                          .nodesize=c(50,100,150),
                          .maxnodes=c(5,10,15,20))
  rf_default <- train(as.factor(Group)~., data=data_boruta,
                      method=customRF,
                      metric="ROC",
                      preProcess = c("center", "scale"),
                      tuneGrid=tunegrid,
                      trControl=control)
  best_rf <- rf_default$bestTune
  #print(rf_default)
  #print("get_model metrix") 
  return(best_rf)
}
tune_model_micro <- tune_model(data_boruta,Group,5)
# ROC was used to select the optimal model using the largest value.
# The final values used for the model were mtry = 4, ntree = 501,
# nodesize = 50 and maxnodes = 20.


get_self_cv <- function(data,Group,best_rf,num){
  set.seed(1234*num)
  fold <- createFolds(y = Group, k=5)
  #return(fold)
  metrics <- matrix(NA,5,13)
  colnames(metrics) <- c("Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision",
                         "Recall","F1","Prevalence","Detection Rate","Detection Prevalence",
                         "Balanced Accuracy","AUC","AUPR")
  for(j in 1:5){
    #print(j)
    fold_test <- data[fold[[j]],]
    fold_train <- data[-fold[[j]],]
    print(table(fold_train$Group))
    Group <- fold_train$Group
    set.seed(1234)
    fold_fit <- randomForest(as.factor(Group)~., data=fold_train,mtry=best_rf$mtry,
                             ntree=best_rf$ntree,importance=TRUE)
    fold_pred <- predict(fold_fit,fold_test)
    #print(fold_pred)
    result <- confusionMatrix(factor(as.vector(fold_pred)),
                              as.factor(fold_test$Group),mode = "everything",
                              positive='Primary_Tumor')
    metrics[j,1:11] <- as.numeric(result$byClass)
    predicted_probs <- predict(fold_fit, fold_test, type = 'prob')
    pred <- prediction(predicted_probs[, 'Solid_Tissue_Normal'], fold_test$Group)
    auc <- performance(pred, 'auc')
    #print(auc@y.values[[1]])
    metrics[j,12] <- auc@y.values[[1]]
    aupr <- performance(pred, 'aucpr')
    metrics[j,13] <- aupr@y.values[[1]]
  }
  # return(metrics)
  print(metrics)
  write.csv(metrics,"result_07_get_self_cv.csv")
  best_index <- which(metrics[,'AUC']==max(metrics[,'AUC']))[1]
  fold_test_best <- data[fold[[best_index]],]
  fold_train_best <- data[-fold[[best_index]],]
  Group <- fold_train_best$Group
  best_model <- randomForest(as.factor(Group)~., data=fold_train_best,mtry=best_rf$mtry,
                             ntree=best_rf$ntree,importance=TRUE)
  result_list <- list("model" = best_model,"metrics" = metrics,"best_fold"=best_index)
  result_list
  # print(result_list)
  # return(result_list)
}
Group <- data_boruta$Group
get_self_cv_micro <- get_self_cv(data_boruta,Group,best_rf,4)

######function-单拎回来#####
set.seed(1234*4)
fold <- createFolds(y = Group, k=5)
#return(fold)
metrics <- matrix(NA,5,13)
colnames(metrics) <- c("Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision",
                       "Recall","F1","Prevalence","Detection Rate","Detection Prevalence",
                       "Balanced Accuracy","AUC","AUPR")
for(j in 1:5){
  #print(j)
  fold_test <- data_boruta[fold[[j]],]
  fold_train <- data_boruta[-fold[[j]],]
  print(table(fold_train$Group))
  Group <- fold_train$Group
  set.seed(1234)
  fold_fit <- randomForest(as.factor(Group)~., data=fold_train,mtry=best_rf$mtry,
                           ntree=best_rf$ntree,importance=TRUE)
  fold_pred <- predict(fold_fit,fold_test)
  #print(fold_pred)
  result <- confusionMatrix(factor(as.vector(fold_pred)),
                            as.factor(fold_test$Group),mode = "everything",
                            positive='Primary_Tumor')
  metrics[j,1:11] <- as.numeric(result$byClass)
  predicted_probs <- predict(fold_fit, fold_test, type = 'prob')
  pred <- prediction(predicted_probs[, 'Solid_Tissue_Normal'], fold_test$Group)
  auc <- performance(pred, 'auc')
  #print(auc@y.values[[1]])
  metrics[j,12] <- auc@y.values[[1]]
  aupr <- performance(pred, 'aucpr')
  metrics[j,13] <- aupr@y.values[[1]]
} 
write.csv(metrics,"result_07_get_self_cv_logcpm_snm.csv")

##虽然4折的auc最高，但5的auc、aupr都很不错
#best_index <- which(metrics[,'AUC']==max(metrics[,'AUC']))[1]
#4

fold_test_best <- data_boruta[fold[[5]],]
fold_train_best <- data_boruta[-fold[[5]],]
Group <- fold_train_best$Group
set.seed(1234)
best_model <- randomForest(as.factor(Group)~., data=fold_train_best,mtry=best_rf$mtry,
                           ntree=best_rf$ntree,importance=TRUE)
result_list <- list("model" = best_model,"metrics" = metrics,"best_fold"=5)
result_list
predicted_probs_best <- predict(best_model,fold_test_best, type = 'prob')
predicted_probs_best2 <- predict(best_model,fold_test_best)

pred_best <- prediction(predicted_probs[, 'Solid_Tissue_Normal'],fold_test_best$Group)

#画图和混淆矩阵
library(pROC)
library(tidyverse)
library(caret)
library(glmnet)
library(patchwork)
library(ROCR)

result_best <- confusionMatrix(factor(as.vector(predicted_probs_best2)),
                          as.factor(fold_test_best$Group),mode = "everything",
                          positive='Primary_Tumor')
result_best
# Reference
# Prediction            Primary_Tumor Solid_Tissue_Normal
# Primary_Tumor                  75                   2
# Solid_Tissue_Normal             0                   4
# 
# Accuracy : 0.9753         
# 95% CI : (0.9136, 0.997)
# No Information Rate : 0.9259         
# P-Value [Acc > NIR] : 0.05536        
# 
# Kappa : 0.7874         
# 
# Mcnemar's Test P-Value : 0.47950        
#                                          
#             Sensitivity : 1.0000         
#             Specificity : 0.6667         
#          Pos Pred Value : 0.9740         
#          Neg Pred Value : 1.0000         
#               Precision : 0.9740         
#                  Recall : 1.0000         
#                      F1 : 0.9868         
#              Prevalence : 0.9259         
#          Detection Rate : 0.9259         
#    Detection Prevalence : 0.9506         
#       Balanced Accuracy : 0.8333         
#                                          
#        'Positive' Class : Primary_Tumor


auc_best <- performance(pred_best,"auc")
perf_auc <- performance(pred_best,"tpr","fpr")

aupr_best <- performance(pred_best, 'aucpr')
perf_aupr <- performance(pred_best,"prec","rec")

plot(perf_auc, 
     col="blue",ylim=c(0,1),
     xlim=c(0,1),
     max.auc.polygon=TRUE,print.auc=T) 
lines(c(0,0),c(1,1),col = "gray", lty = 4 )
text(0.8,0.6, labels = paste0("AUC = ",round(metrics[5,12],3)))

plot(perf_aupr, 
     col="blue",ylim=c(0,1),
     xlim=c(0,1),
     max.auc.polygon=TRUE,print.auc=T) 
lines(c(0,0),c(1,1),col = "gray", lty = 4 )
text(0.6,0.6, labels = paste0("AUPR = ",round(metrics[5,13],3)))

all(rownames(otu_table_logcpm_snmdata)==rownames(sample_table_fil))
primay_name <- read.csv("/home/yuekaile/keller/stad_rna/sample-name.csv",
                        header = T,row.names = 1,check.names = F)
primay_name$x <- gsub(".","-",primay_name$x,fixed = T)
otu_table_logcpm_snmdata2 <- otu_table_logcpm_snmdata[primay_name$x,]

write.csv(otu_table_logcpm_snmdata,file = "otu_table_logcpm_snmdata_407.csv")
write.csv(otu_sig,file = "otu_table_logcpm_snmdata_sigspecies_407.csv")
write.csv(otu_table_logcpm_snmdata2,file = "otu_table_logcpm_snmdata_338.csv")

data_boruta2 <- data_boruta
data_boruta2 <- data_boruta2[,rownames(feature_importance)]
boxplot(data_boruta2,col = "red")

boruta.variable.imp.use2 <- boruta.variable.imp.use[,rownames(feature_importance)]
boxplot(boruta.variable.imp.use2,col="red",widths=0.8)
