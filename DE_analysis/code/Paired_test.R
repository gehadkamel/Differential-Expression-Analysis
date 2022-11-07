#Paired test function 
# METHODE 1 : Hypothesis testing if data are paired 
#If paired difference follow normal distribution, use t-test
#Ho : GE in healthy tissues = GE in cancer tissues
#H1 : GE in healthy tissues != GE in cancer tissues
t_test <- function (data.healthy, data.cancer)
{
  #data.healthy <- GE.healthy_subset
  #data.cancer <- GE.cancer_subset
  paired_diff <- data.healthy - data.cancer
  GN_data <- row.names(data.healthy)
  #Check for paired difference normality
  paired_p_value = vector(mode = "double", length = nrow(paired_diff))
  for (i in 1:nrow(paired_diff))
  {
    
    paired_p_value[i] <- shapiro.test(as.matrix(paired_diff[i,]))$p.value
    
  }
  test_p_value = vector(mode = "double", length = nrow(paired_diff))
  #performing paired test
  for (i in 1:nrow(data.healthy))
  {
    if (paired_p_value[i] > 0.05)
    {
      #data follow normal distribution, perform t-test
      test_p_value[i] <- t.test(as.matrix(as.numeric(data.healthy[i,])),as.matrix(as.numeric(data.cancer[i,])), paired = TRUE, alternative = "two.sided" )$p.value
    }
    else
    {
      #data do not follow normal distribution, peroform wilcox test 
      test_p_value[i] <- wilcox.test(as.matrix(as.numeric(data.healthy[i,])), as.matrix(as.numeric(data.cancer[i,])), paired = TRUE , alternative = "two.sided")$p.value
    }
    
  }
  paired_test_p_value = data.frame(GN_data , test_p_value)
  
  return(paired_test_p_value)
  
}


#Fold change function
fold_change <- function(healthy.data, cancer.data)
{
  #healthy.data <- GE.healthy_lusc
  #cancer.data <- GE.cancer_lusc
  GN_fc <- row.names(healthy.data)
  #Getting the mean of each  gene across the  samples
  healthy.data <- apply(healthy.data, 1, mean)
  cancer.data <- apply(cancer.data,1 , mean)
  
  #getting fold change between healthy and cancer data
  fc_log2 <- log2(healthy.data) - log2(cancer.data)
  
  #Defining the threshold to identify DEGs
  threshold <- log2(1.5)
  #getting absolute values to compare fold change 
  fc_log2_abs <- abs(fc_log2)
  fc_log2_abs <- as.data.frame(fc_log2_abs)  
  #Define colnames for empty data frames
  
  x <- c("Fold change", "Gene")
  fold_DEGs <- data.frame(matrix(ncol = 2, nrow = 0))
  fold_not_DEGs <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(fold_DEGs) <- x
  colnames(fold_not_DEGs) <- x
  #Looping to differentiate between DEGs and
  for (i in 1:length(fc_log2))
 {
  #if the fold change > threshold then, it is DEG
  if (fc_log2_abs[i,] >= threshold)
  {
    fold_DEGs[i,1] <- fc_log2_abs[i,]
    fold_DEGs[i,2] <- GN_fc[i]
  }
  #If less then the threshold , then not DEG
  else
  {
    fold_not_DEGs[i,1] <- fc_log2_abs[i,]
    fold_not_DEGs[i,2] <- GN_fc[i]
  }
 }
  fold_DEGs <- na.omit(fold_DEGs)
  fold_not_DEGs <- na.omit(fold_not_DEGs)
  return(list(fold_DEGs, fold_not_DEGs,fc_log2))
  

}


###########################################################
run_enhanced_volcano <- function(Data_Healthy, Data_Cancer ,res,cancer_type){
  
  
  #volcano plot Data Preparation
  #pvalue <- run_paired(Data_Healthy, Data_Cancer, cancer_type)
  #log2FoldChange <- run_fc(Data_Healthy, Data_Cancer,cancer_type)
  #res <- data.frame( log2FoldChange,pvalue)
  
  #print(paste("Starting Volcano Plot for",cancer_type))
  
  # The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
  
  out <- EnhancedVolcano(res,
                         lab = rownames(res),
                         x = 'FC',
                         y = 'rawPvalue',
                         pCutoff = 0.05,   #alpha - on Y_axis
                         FCcutoff = 0.6,    #Threshold - on X_axis
  )
  plot(out)
  
  print("Ending  Volcano Plot .")
  
  #Extraction  Names of DEGS Using volcano plot
  Volcano_DEGS <- out$layers[[4]]$data$lab
  
  #Preparation DEGS For GSEA Software
  Volcano_DEGS_DF_Healthy <- Data_Healthy[rownames(Data_Healthy) %in% Volcano_DEGS, ] 
  Volcano_DEGS_DF_Cancer <- Data_Cancer[rownames(Data_Cancer) %in% Volcano_DEGS, ]
  DEGS_ChangeN_Cancer <- sub("T", "H", colnames(Volcano_DEGS_DF_Cancer)) 
  colnames(Volcano_DEGS_DF_Cancer) <-DEGS_ChangeN_Cancer
  Volcano_DEGS_DH <-data.frame( Volcano_DEGS_DF_Cancer,Volcano_DEGS_DF_Healthy )
  write.table(Volcano_DEGS_DH, paste(cancer_type,"Volcano_DEGs.txt",sep="_"),row.names = TRUE,sep = "\t")
  
}