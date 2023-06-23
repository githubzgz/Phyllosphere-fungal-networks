######################################################
######Contribution of MEN modules to saponin(H2)######
######################################################
#Ctrl+Shift+O:catalogue of this script
#R version 4.2.2

#Load required packages
library(dplyr)
library(stringr)
library(Matrix)
library(igraph)
library(reshape2)
library(tidyverse)
library(mlr3)
library(mlr3extralearners)
library(mlr3verse)
library(vegan)

#1.Calculating module eigengenes-----------
#Can also use the moduleeigengenes function in WGCNA package
#(1)Extracting modules with number of nodes higher than 5
module_info <- rbind(aggregate(num~greedy.module,data=f_gin_le_ig_node%>%mutate(.,num=1),FUN=sum) %>% .[.$num>=5,]%>%mutate(.,plant='gin',compart='_le'),
                     aggregate(num~greedy.module,data=f_gin_lp_ig_node%>%mutate(.,num=1),FUN=sum) %>% .[.$num>=5,]%>%mutate(.,plant='gin',compart='_lp'),
                     
                     aggregate(num~greedy.module,data=f_qui_le_ig_node%>%mutate(.,num=1),FUN=sum) %>% .[.$num>=5,]%>%mutate(.,plant='qui',compart='_le'),
                     aggregate(num~greedy.module,data=f_qui_lp_ig_node%>%mutate(.,num=1),FUN=sum) %>% .[.$num>=5,]%>%mutate(.,plant='qui',compart='_lp'),
                     
                     aggregate(num~greedy.module,data=f_not_le_ig_node%>%mutate(.,num=1),FUN=sum) %>% .[.$num>=5,]%>%mutate(.,plant='not',compart='_le'),
                     aggregate(num~greedy.module,data=f_not_lp_ig_node%>%mutate(.,num=1),FUN=sum) %>% .[.$num>=5,]%>%mutate(.,plant='not',compart='_lp'))

#(2)Calculating eigengenes
my_module_eigengenes <- function(module_info,plant,compart,ig_node,comm){
  require(vegan)
  n = module_info[module_info$plant==plant & module_info$compart==compart,'greedy.module'] %>% as.vector %>% as.numeric
  result <- data.frame(matrix(NA,nrow=27,ncol = length(n)))
  names(result) = paste(paste(plant,compart,sep = ''),n,sep='_module')
  for(i in 1:length(n)){
    module = n[i]
    node = ig_node[ig_node$greedy.module==module,'Row.names']
    module_comm <- comm[rownames(comm)%in% node,1:27]
    module_pca <- rda(module_comm%>%t,scale=F)
    module_pca <- module_pca$CA$u[,1]%>%as.vector%>%as.numeric
    result[,i]=module_pca
    
  }
  return(result)
}

#2.Machine learning framework-------------
#a.PCA of saponin contents(Using PG as an example)
gin_le_sap_pc1 = rda(gin_le_saponin[,1:8],scale=F)$CA$u[,1]#PC1:59.07%

#b.Machine learning evaluation
my_ml <- function(ml_name,
                  bs_factor,
                  saponin_pc1,
                  le_module,lp_module){
  #ml_name:key value of machine learners;
  #bs_factor:edaphic factors with no strong collinearity;
  #saponin:saponin contents in leaf;
  #le_module:module eigengenes of le(Only candidate fungal modules);
  #lp_module:module eigengenes of lp(Only candidate fungal modules).
  require(mlr3extralearners)
  require(mlr3verse)
  
  #(1)Generating tasks
  env_data = bs_factor %>% mutate(.,sap_pc1=saponin_pc1$pc1)
  env_le_data = cbind(env_data,le_module)
  env_lp_data = cbind(env_data,lp_module)
  env_all_data = cbind(env_data,le_module,lp_module)
  
  env_task = as_task_regr(env_data,target='sap_pc1',id='env')
  env_le_task = as_task_regr(env_le_data,target='sap_pc1',id='env_le')
  env_lp_task = as_task_regr(env_lp_data,target='sap_pc1',id='env_lp')
  env_all_task = as_task_regr(env_all_data,target='sap_pc1',id='all')
  
  #(2)Regression using different learners
  result = data.frame()
  for(i in 1:length(ml_name)){
    seed=sample(1:1e06,1)#Selecting random seed
    #a.Creating learners
    my_lrn = lrn(ml_name[i])
    #b.Creating resampling strategy
    my_resample = rsmp('cv',
                       folds=5)
    #c.Resampling analysis
    set.seed(seed)
    env_resample = resample(env_task,my_lrn,my_resample)
    env_le_resample = resample(env_le_task,my_lrn,my_resample)
    env_lp_resample = resample(env_lp_task,my_lrn,my_resample)
    env_all_resample = resample(env_all_task,my_lrn,my_resample)
    #d.Creating performance evaluation methods
    my_measure = msrs(c('regr.mse',#Mean Squared Error
                        'regr.srho',#Spearman's rho))
    #e.Performance assessment
    env_resample_measure = env_resample$score(my_measure)
    env_le_resample_measure = env_le_resample$score(my_measure)
    env_lp_resample_measure = env_lp_resample$score(my_measure)
    env_all_resample_measure = env_all_resample$score(my_measure)
    #f.Combining results
    all_measure = rbind(env_resample_measure,env_le_resample_measure,env_lp_resample_measure,env_all_resample_measure)%>%mutate(.,seed=seed)
    result <- rbind(result,all_measure)
  }
  return(result)
}

#b.Appointing learners
ml_name = c('regr.IBk',#K-nearest neighbour
            'regr.randomForest',#RandomForest
)
#c.An example using PG
gin_ml_pc1_result <- my_ml(ml_name,gin_bs_factor_no_collinear,
                           gin_le_sap_pc1,
                           f_gin_ig_module_eigengene,
                           f_gin_ig_module_eigengene)


#3.SEM construction----------------
gin_sem_data = cbind(
  gin_bs_factor_no_collinear,#Variables with strong collinearity have been removed
  f_gin_ig_module_eigengene[,c('gin_le_module1','gin_le_module2','gin_lp_module1','gin_lp_module2','gin_lp_module3')])%>%
  mutate(.,le_sum=gin_le_saponin$sum)

gin_sap_model <- '
#Regression
gin_le_module2 ~ OC+AK+NIN+MC+Ca
gin_lp_module3 ~ OC+AK+NIN+MC+Ca
le_sum ~ gin_le_module2+gin_lp_module3+OC+AK+MC+NIN

#Correlation
OC ~~ AK+NIN+MC+Ca
AK ~~ NIN+MC+Ca
NIN ~~ MC+Ca
MC ~~ Ca

gin_le_module2~~gin_lp_module3

'
gin_sap_fit <- sem(gin_sap_model,data=gin_sem_data%>%scale)
fitMeasures(gin_sap_fit,c("chisq","df","pvalue","cfi","rmsea"))
summary(gin_sap_fit,standardized=T)
