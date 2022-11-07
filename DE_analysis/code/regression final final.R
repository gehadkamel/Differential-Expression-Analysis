library(glmnet)
library(car)
setwd ("E:/bioinformatics/statistics/projects/GEHAD") #open data files

#######################################
#function feature selection
run_Feature_Selection <- function(Data.CNV , One.DEGS) {
  fit.cv <- cv.glmnet(Data.CNV, One.DEGS, family="gaussian", alpha=1, standardize=FALSE, nfolds=5)
  # Getting the value of lowest value of lambda.
  lambda <- fit.cv$lambda.min  
  # Here we compute the variables regression coefficeints after being penalized (set to zero).
  model <- glmnet(Data.CNV, One.DEGS, family="gaussian", alpha=1, lambda=lambda, standardize=FALSE)
  # exclude intercept (the first value).
  coef.fit <- coef(model, s=lambda)[2:(ncol(Data.CNV)+1)]
  # Keep features with non-zero regression coefficients.
  features.in <- which(abs(coef.fit) > 0)
  # The variables matrix after removing the penalized features.
  FS.CNV = Data.CNV[,features.in]
  return(FS.CNV)
}
############################################
# regression for lusc 
# open lusc data files  
top_5_lusc=read.table("Top_5_lusc_DEGs.csv", sep=",", header = TRUE , row.names= 1)
cnv_data_o=read.table("CNV_lusc.csv", sep=",", header = TRUE , row.names= 1)
CNV_lusc_lm=read.table("CNV_lusc_lm.csv", sep=",", header = TRUE , row.names= 1)

top_5_lusc= t(top_5_lusc)
CNV_lusc_lm= t(CNV_lusc_lm)
t = row.names(CNV_lusc_lm)
t2=row.names(top_5_lusc)
top5_degs_sample= subset( top_5_lusc, t2 %in% t) #get the 13 common sample of top 5 DEGS
top5_degs_sample=t(top5_degs_sample)

top5_degs_sample=as.matrix(top5_degs_sample)
cnv_data=as.matrix((cnv_data_o))

for ( i in 1:nrow (top5_degs_sample)){
  
  if (ncol(cnv_data) > nrow(cnv_data)){
    
    #Feature Selection --Lung
    cnv_data_1 =  run_Feature_Selection(cnv_data , top5_degs_sample[i,])
    print(row.names(top5_degs_sample)[i])
    # print(top5_degs_sample[i,])
    # print(cnv_data_1)
    cnv_data =cnv_data_1
  }
  # deg_1= top5_degs_sample[i,] 
  # model= lm(deg_1 ~ cnv_data) #linear regression
  model= lm( top5_degs_sample[i,]~ cnv_data) #linear regression
  print(summary(model))
  capture.output(summary(model), file="linear regression for lusc.txt", append =T )
  require(car)
  avPlots(model)
  
  cnv_data=as.matrix((cnv_data_o))
}

######################################
#  Regression for Kirc 

# open Kirc data files  
top_5_kirc =read.table("Top_5_kirc_DEGs.csv", sep=",", header = TRUE , row.names= 1)
cnv_Kirc_data_o =read.table("CNV_kirc.csv", sep=",", header = TRUE , row.names= 1)

top5_degs_Kirc_sample=as.matrix(top_5_kirc)
cnv_Kirc_data_o=as.matrix(cnv_Kirc_data_o)

for ( i in 1:nrow (top5_degs_Kirc_sample)){
  
  if (ncol(cnv_Kirc_data_o) > nrow(cnv_Kirc_data_o)){
    
    #Feature Selection --Lung
    cnv_kirc_data_1 =  run_Feature_Selection(cnv_Kirc_data_o , top5_degs_Kirc_sample[i,])
    print(row.names(top5_degs_Kirc_sample)[i])
    cnv_kirc_data =cnv_kirc_data_1
  }
  model= lm( top5_degs_Kirc_sample [i,] ~ cnv_kirc_data) #linear regression
  print(summary(model))
  capture.output(summary(model), file="linear regression for kirc_f.txt", append =T )
  require(car)
  avPlots(model)
  
}

 