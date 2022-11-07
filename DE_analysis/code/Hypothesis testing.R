#Requirment 1
library(stringr)

#Method1 : Hypothesis testing
#Ho: GE in healthy tissues = GE in cancer tissues
#H1: GE in healthy tissues != GE in cancer tissues
setwd("D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project 1")
#Paths for lung cancer
healthy.path_lusc = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project_Data/lusc-rsem-fpkm-tcga_paired.txt"
cancer.path_lusc = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project_Data/lusc-rsem-fpkm-tcga-t_paired.txt"
lusc_CNV.path = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project_Data/lusc_CNV_core.txt"
#Paths for kidney cancer data

healthy.path_kirc = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project_Data/kirc-rsem-fpkm-tcga_paired.txt"



cancer.path_kirc = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project_Data/kirc-rsem-fpkm-tcga-t_paired.txt"
kirck_CNV.path = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project_Data/kirc_CNV_core.txt"

#Sources for functions 
source("Data prep.R")
source("Paired_test.R")
source("indep_testt.R")

#Data preparation for LUSC and kirc
Data_lusc <- prepare_GE_data (healthy.path_lusc, cancer.path_lusc, lusc_CNV.path)
Data_kirc <- prepare_GE_data (healthy.path_kirc, cancer.path_kirc, kirck_CNV.path)

#Lusc data
GE.healthy_lusc  = as.data.frame(Data_lusc[1])
GE.cancer_lusc = as.data.frame(Data_lusc[2])
CNV_lusc = as.data.frame(Data_lusc[3])
write.csv(CNV_lusc, file = "CNV_lusc.csv")
CNV_lusc_lm = as.data.frame(Data_lusc[4])
write.csv(CNV_lusc_lm, file = "CNV_lusc_lm.csv")
#Paired test for lusc
paired_t_lusc = t_test(GE.healthy_lusc, GE.cancer_lusc)
#indep test for lusc``
indep_t_lusc = indep_test(GE.healthy_lusc , GE.cancer_lusc)
#Pvalue correction
#p.ind_value_adjusted  = p.adjust(indep_t_lusc$indep_test_p_value, method = "fdr")
#p.paired_value_adjusted = p.adjust(paired_t_lusc$test_p_value, method = 'fdr')

#paired_t_lusc$test_p_value <- p.paired_value_adjusted

#indep_t_lusc$indep_test_p_value <- p.ind_value_adjusted

#Identify DEGs from paired test
#Getting genes with p values less than 0.05
DEGs_lusc_ind <- which(paired_t_lusc$test_p_value < 0.05)
DEGs_lusc_healthy <- GE.healthy_lusc[DEGs_lusc_ind,]

DEGs_lusc_cancer <- GE.cancer_lusc[DEGs_lusc_ind,]
GN_paired <- row.names(DEGs_lusc_healthy)
#Identifying DEGs from independent test and plotting p_values 
DEGs_lusc_ind_ind <- which(indep_t_lusc$indep_test_p_value < 0.05 )
DEGs_ind_lusc_healthy <- GE.healthy_lusc[DEGs_lusc_ind_ind, ]
DEGs_ind_lusc_cancer <- GE.cancer_lusc[DEGs_lusc_ind_ind,]
GN_indep <- row.names(DEGs_ind_lusc_healthy)

path_lusc_his <- "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project 1/LUSC_pvlaues.png"

png(path_lusc_his, width  = 800, height = 1000, res = 200)
hist(indep_t_lusc$indep_test_p_value, ylab = "Frequence", xlab = "P-values")
dev.off()


#Identifying common Genes from indep and paired tests
common_lusc <- intersect(row.names(DEGs_lusc_cancer), row.names(DEGs_lusc_cancer))
diff_ind <- setdiff(GN_indep, common_lusc)
diff_paired <- setdiff(GN_paired, common_lusc)


#identify top 5 DEGs from paired test
sorted_p.value = sort(paired_t_lusc$test_p_value, index.return = T, decreasing = F)  #$ix  
#extract the top 5 genes with least p.value 
top_5_p.value = as.data.frame(sorted_p.value)[,2][1:5]
# return the expression values of the most Deferentially expressed genes from the original data 
top_5 = GE.cancer_lusc[top_5_p.value,]
top_5_lusc_lm <- as.matrix(CNV_lusc_lm[top_5_p.value,])
path_name = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project 1"
write.csv(top_5_lusc_lm, file = "Top_5_lusc_DEGs.csv", row.name = T )



#######
#Plotting pie chart identifying DEGs, non DEGs, and not_expressed
names_lusc = c("DEGs", "not_DEGs", "not_expressed")
non_sig_lusc <- nrow(GE.healthy_lusc) - nrow(DEGs_ind_lusc_healthy)
non_DEGs_lusc <- GE.healthy_lusc[-DEGs_lusc_ind,]
numbers_lusc = c(nrow(DEGs_ind_lusc_cancer), nrow(non_DEGs_lusc), non_sig_lusc)
percentage_lusc = round(numbers_lusc*100/sum(numbers_lusc),1)
percen = paste(percentage_lusc,"%", sep = "")
labels_lusc = paste(names_lusc,":", percen, sep = "")
colors = c("blue", "red", 'black')
path_lusc <- "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project 1/Lusc.png"
png(path_lusc, width = 1000, height = 600, res = 200)

pie(numbers_lusc, labels = labels_lusc, col = colors , main = "LUSC")

dev.off()
###########################################################
#Method 2: Fold change
fold_change_lusc <- fold_change(GE.healthy_lusc, GE.cancer_lusc)
DEGs_fold_lusc <- as.data.frame(fold_change_lusc[1])
not_DEGs_fold_lusc <- as.data.frame(fold_change_lusc[2])
fc_log_lusc <- as.data.frame(fold_change_lusc[3])
colnames(fc_log_lusc) <- "FC"
write.csv(DEGs_fold_lusc, file = "Fold_DEGs_lusc.csv")
write.csv (not_DEGs_fold_lusc, file = "Fold_not_DEGs_lusc.csv")
#Plotting fold change in lusc
path_fc_lusc <- "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project 1/Lusc_fc.png"
png(path_fc_lusc, width = 800, height = 800, res = 200)
hist(fc_log_lusc$FC, xlab = 'Gene fold change',ylim = c(0,6000),main = "Fold change in lusc")
dev.off()
############################################################
#Method3 : Volcano plotting
#Compare results between fold change and paired testing 

#P values from paired testing
p.paired_value_adjusted
#fold change from FC
fc_log_lusc
#creating data frame containing rawpvalues and fc
result <- data.frame(as.numeric(fc_log_lusc[,1]), paired_t_lusc$test_p_value)
colnames(result) <- c("FC", "rawPvalue")
rownames(result) <- rownames(GE.healthy_lusc)
write.csv(result, file = "fc_paired_results.csv")
#Volcano plotting

run_enhanced_volcano(GE.healthy_lusc, GE.cancer_lusc, result, "LUSC")


########################################################################################
##########################################################################################

#KIRC data 
healthy.path_kirc = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project_Data/kirc-rsem-fpkm-tcga_paired.txt"
cancer.path_kirc = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project_Data/kirc-rsem-fpkm-tcga-t_paired.txt"
cnv.path_kirc  = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project_Data/kirc_CNV_core.txt"

#Data preparation for KIRC

Data_kirc <- prepare_GE_data(healthy.path_kirc, cancer.path_kirc, cnv.path_kirc)

#KIRC data
GE.healthy_kirc  = as.data.frame(Data_kirc[1])
GE.cancer_kirc = as.data.frame(Data_kirc[2])
CNV_kirc = as.data.frame(Data_kirc[3])
write.csv(CNV_kirc, file = "CNV_kirc.csv")
CNV_kirc_lm = as.data.frame(Data_kirc[4])
write.csv(CNV_kirc_lm, file = "CNV_lusc_lm.csv")
#Paired test for KIRC
paired_t_kirc = t_test(GE.healthy_kirc, GE.cancer_kirc)
#indep test for KIRC
indep_t_kirc = indep_test(GE.healthy_kirc , GE.cancer_kirc)
#Pvalue correction
#p.ind_value_adjusted_kirc  = p.adjust(indep_t_kirc$indep_test_p_value, method = "fdr")
#p.paired_value_adjusted_kirc = p.adjust(paired_t_kirc$test_p_value, method = 'fdr')

#paired_t_kirc$test_p_value <- p.paired_value_adjusted_kirc

#indep_t_lusc$indep_test_p_value <- p.ind_value_adjusted_kirc

#Identify DEGs from paired test
#Getting genes with p values less than 0.05
DEGs_kirc_ind <- which(paired_t_kirc$test_p_value < 0.05)
DEGs_kirc_healthy <- GE.healthy_kirc[DEGs_kirc_ind,]

DEGs_kirc_cancer <- GE.cancer_kirc[DEGs_kirc_ind,]
GN_paired_kirc <- row.names(DEGs_kirc_healthy)
#Identifying DEGs from independent test
DEGs_kirc_ind_ind <- which(indep_t_kirc$indep_test_p_value < 0.05 )
DEGs_ind_kirc_healthy <- GE.healthy_kirc[DEGs_kirc_ind_ind, ]
DEGs_ind_kirc_cancer <- GE.cancer_kirc[DEGs_kirc_ind_ind,]
GN_indep_kirc <- row.names(DEGs_ind_kirc_healthy)

path_kirc_his <- "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project 1/KIRC_pvlaues.png"

png(path_kirc_his, width  = 800, height = 1000, res = 200)
hist(indep_t_kirc$indep_test_p_value, ylab = "Frequence", xlab = "P-values")
dev.off()
#Identifying common Genes from indep and paired tests
common_kirc <- intersect(row.names(DEGs_kirc_cancer), row.names(DEGs_kirc_cancer))
diff_ind_kirc <- setdiff(GN_indep_kirc, common_kirc)
diff_paired_kirc <- setdiff(GN_paired_kirc, common_kirc)
#identify top 5 DEGs from paired test
sorted_p.value_kirc = sort(paired_t_kirc$test_p_value, index.return = T, decreasing = F)  #$ix  
#extract the top 5 genes with least p.value 
top_5_p.value_kirc = as.data.frame(sorted_p.value_kirc)[,2][1:5]
# return the expression values of the most Deferentially expressed genes from the original data 
top_5_kirc = GE.cancer_kirc[top_5_p.value_kirc,]
top_5_kirc_lm <- as.matrix(CNV_kirc_lm[top_5_p.value_kirc,])
path_name_kirc = "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project 1"
write.csv(top_5_kirc_lm, file = "Top_5_kirc_DEGs.csv", row.name = T )




###########################################################
#Method 2: Fold change
fold_change_kirc <- fold_change(GE.healthy_kirc, GE.cancer_kirc)
DEGs_fold_kirc <- as.data.frame(fold_change_kirc[1])
not_DEGs_fold_kirc <- as.data.frame(fold_change_kirc[2])
fc_log_kirc <- as.data.frame(fold_change_kirc[3])
colnames(fc_log_kirc) <- "FC"
write.csv(DEGs_fold_kirc, file = "Fold_DEGs_kirc.csv")
write.csv (not_DEGs_fold_kirc, file = "Fold_not_DEGs_kirc.csv")

#plotting fold change 
path.kirc_fc <- "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project 1/KIRC_fc.png"

png(path.kirc_fc , width = 800, height = 1000, res = 200)
hist(fc_log_kirc$c.METTL21B....0.335404617282969..TXNDC16....0.0030790929951543.., xlab = "Fold change values", main = "KIRC fold change of DEGs")

dev.off()
####################################################################
#Plotting pie chart identifying DEGs, non DEGs, and not_expressed
names_kirc = c("DEGs", "not_DEGs", "not_expressed")
non_sig_kirc <- nrow(GE.healthy_kirc) - nrow(DEGs_kirc_ind)
non_DEGs_kirc <- GE.healthy_kirc[-DEGs_kirc_ind,]
numbers_kirc = c(nrow(DEGs_ind_kirc_cancer), nrow(non_DEGs_kirc), non_sig_kirc)
percentage_kirc = round(numbers_kirc*100/sum(numbers_kirc),1)
percen_kirc = paste(percentage_kirc,"%", sep = "")
labels_kirc = paste(names_kirc,":", percen_kirc, sep = "")
colors = c("blue", "red", 'black')
path_KIRC <- "D:/Bioinformatics/Nile University/Statistical Analysis/Course project/Project 1/KIRC.png"
png(path_KIRC, width = 1000, height = 600, res = 200)

pie(numbers_kirc, labels = labels_kirc, col = colors , main = "KIRC")

dev.off()



#Volcano plotting 
result_kirc <- data.frame(as.numeric(fc_log_kirc[,1]), paired_t_kirc$test_p_value)
colnames(result_kirc) <- c("FC", "rawPvalue")
rownames(result_kirc) <- rownames(GE.healthy_kirc)
write.csv(result_kirc, file = "fc_paired_kirc_results.csv")


run_enhanced_volcano(GE.healthy_kirc, GE.cancer_kirc, result_kirc, "KIRC")



