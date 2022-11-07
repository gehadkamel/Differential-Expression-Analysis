#Regression and feature selection
#Predictor: CNV , Response: GE of specific genes across tissues


run_reg <- function(top_5_DEGs , data.cnv)
{
  
  top_5_DEGs <- as.matrix(top_5_DEGs)
  data.cnv <- as.matrix(data.cnv)
  
  for (i in 1:nrow(top_5_DEGs))
  {
    if(ncol(data.cnv) > nrow(data.cnv))
    {
      #Run feature selection if number of variables(predictors) more than data points(samples)
      data.cnv <- feature_selection(data.cnv, top_5_DEGs[i,])
      
    }
    
    DEGs.1 <- top_5_DEGs[i,]
    #performing linear modeling whether feature selection occured or number of variables were less than number of 
    #responses
    model <- lm(DEGs.1 ~ data.cnv)
    print(summary(model))
    #Getting results into file 
    capture.output(summary(model),file = "Linear_regression_Coef.txt", append = T)
    require(car)
    avPlots(model)
    
    
  }
  
  
}

feature_selection <- function (Data.CNV, DEGs_1)
{
  #Convert data into matrix
  Data.CNV <- as.matrix(Data.CNV)
  #Run cross validation to obtain min lambda values
  fit.cv <- cv.glmnet(Data.CNV, DEGs_1, family = "gaussian", alpha = 1, standardize = FALSE, nfolds = 5)
  #obtain the minimum lambda value
  lambda <- fit.cv$lambda.min
  #Compute value of the variable regression coefficient after being penalized and removing zero
  model <- glmnet(Data.CNV, DEGs_1, family = "gaussian", alpha = 1 ,standardize = FALSE)
  #Exclude intercept %first value
  coef.fit <- coef(model, s = lambda)[2:ncol(Data.CNV+1)]
  #Getting the indices of the significant features(predictors)
  features.in <- which(abs(coef.fit) > 0)
  #Obtaining matrix after removing non significant features
  FS.CNV <- Data.CNV[,features.in]
  
  return(FS.CNV)

  
}
