#Method 1: hypothesis testing
#Data are independent 
#Ho : GE in healthy tissues = GE in cancer tissues
#H1 : GE in healthy tissues != GE in cancer tissues
indep_test <- function(data.healthy, data.cancer)
{
  
  Genes_data <- row.names(data.healthy)
  indep_test_p_value = vector(mode = "double", length = nrow(data.healthy))
  for (i in 1:nrow(data.healthy))
  {
    #Convert data into  numeric 
    #check for normality for each subset of data
    if (shapiro.test(as.matrix(as.numeric(data.healthy[i,])))$p.value > 0.05 && shapiro.test(as.matrix(as.numeric(data.cancer[i,])))$p.value > 0.05)
    
    {
      #if both follow normal distribution then we perform variance check
      if(var.test(as.matrix(as.numeric(data.healthy[i,])),as.matrix(as.numeric(data.cancer[i,])))$p.value > 0.05)
      {
        indep_test_p_value[i] <- t.test(as.matrix(as.numeric(data.healthy[i,]),as.matrix(as.numeric(data.cancer[i,])), paired = FALSE, var.equal = TRUE , alternative = "two.sided"))$p.value
      }
      #if data variances are not equall
      else
      {
        indep_test_p_value[i] <- t.test(as.matrix(as.numeric(data.healthy[i,])),as.matrix(as.numeric(data.cancer[i,])), paired = FALSE, var.equal = FALSE , alternative = "two.sided")$p.value
      
      }
    }
    #if one of the data does not follow normal distribution, perform wilcoxon rank sum test  
    else
    {
      indep_test_p_value[i] <- wilcox.test(as.matrix(as.numeric(data.healthy[i,])), as.matrix(as.numeric(data.cancer[i,]), alternative = "two.sided"))$p.value
      
    }
      
    
  }
  return (indep_test_p_value <- data.frame(Genes_data, indep_test_p_value))
  
}

