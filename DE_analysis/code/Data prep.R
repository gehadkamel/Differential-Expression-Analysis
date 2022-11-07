#Data preparation

prepare_GE_data <- function(path.healthy, path.cancer,path.CNV)
{

  #Read GE data 
  GE.healthy <- read.delim(path.healthy, header = T, sep = "\t")
  GE.cancer <- read.delim(path.cancer , header =T , sep = "\t")
  #defining rownames as gene names 
  row.names(GE.healthy) <- GE.healthy [,1]
  row.names(GE.cancer)  <- GE.cancer [,1]
  
  #removing extraGene names column and Gene IDs 
  GE.healthy <- GE.healthy[, -1:-2]
  GE.cancer <- GE.cancer [, -1:-2]

  #Exclude genes that 50% of  values == 0 
  filter.GE.healthy <- GE.healthy [which(apply(GE.healthy, 1, median)!=0), ] 
  filter.GE.cancer <- GE.cancer[which(apply(GE.cancer, 1 , median) != 0), ]
  
  #Getting gene names after filteration
  
  healthy_gene.names <- rownames(filter.GE.healthy)
  cancer_gene.names <- rownames(filter.GE.cancer)
  common <- intersect(healthy_gene.names, cancer_gene.names)
  
  #Extraction of common gene names in both data
  GE.healthy_subset <- filter.GE.healthy[rownames(filter.GE.healthy) %in% common,]
  GE.cancer_subset <- filter.GE.cancer[rownames(filter.GE.cancer) %in% common,] 
  #Read CNV data
  CNV_data <- read.delim(path.CNV, header =TRUE , sep = "\t")
  
  #Filter data with zero values and Getting sample names as present in GE data
  filter.CNV <- CNV_data[, which(apply(CNV_data, 2, median)!=0)]
  row.names(filter.CNV) <- CNV_data[,1]
  filter.CNV <- filter.CNV[,-1]
  CNV.F_rowname <- row.names(filter.CNV)
  CNV.F.replace_rowname<- str_replace_all(CNV.F_rowname,"-",".")
  row.names(filter.CNV) <- CNV.F.replace_rowname
  CNV.F_Arowname <- row.names(filter.CNV)
  
  #Getting common samples between CNV and GE data
  cancer_samples <- colnames(GE.cancer_subset)
  Common_CNV <- intersect(cancer_samples, CNV.F_Arowname)
  
  data.CNV <- filter.CNV[rownames(filter.CNV) %in% Common_CNV, ]
  data_cancer_lm <- GE.cancer_subset[ ,colnames(GE.cancer_subset) %in% Common_CNV]
  
  return(list(GE.healthy_subset,GE.cancer_subset,data.CNV, data_cancer_lm))
  
}
