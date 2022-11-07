library(ggplot2)
library(dplyr)
library(ggrepel)

LUSC_FC_P= read.csv("D:\\Bioinformatics Diploma\\Statistical Analysis Project\\Project Code\\fc_paired_results.csv",
                    header = TRUE)
names(LUSC_FC_P)[1] <- 'HUGO'

KIRC_FC_P= read.csv("D:\\Bioinformatics Diploma\\Statistical Analysis Project\\Project Code\\fc_paired_kirc_results.csv", 
                      header = TRUE)
names(KIRC_FC_P)[1] <- 'HUGO'

#Creating Volcano Plot for LUSC 

#Adding a new column to label values to make it easier for plotting
LUSC_FC_P$DifferentiallyExpressed <-"NO Difference"  

#Setting the threshold for p-value and Log fold change 
LUSC_FC_P$DifferentiallyExpressed[LUSC_FC_P$FC > 0.585 & LUSC_FC_P$rawPvalue < 0.05] <- "Upregulated"
LUSC_FC_P$DifferentiallyExpressed[LUSC_FC_P$FC < -0.585 & LUSC_FC_P$rawPvalue < 0.05] <- "Downregulated"

#Using ggplot for plotting and ggrepel for labeling
ggplot(data=LUSC_FC_P, aes(x=FC, y=-log10(rawPvalue), col=DifferentiallyExpressed, label=HUGO )) +
  geom_point()+                #Creating Scatter plot
  theme_minimal()+             #Theme with no background annotaions 
  geom_label_repel()+           #Used to avoid overlapping of labels 
  geom_text_repel()+
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-0.585, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")+
  coord_cartesian(xlim =c(-10, 10), ylim = c(0, 20))+ 
  ggtitle("LUSC Volcano plot")    #Labeling the figure 


#############################################################################################################

#Creating Volcano Plot for KIRC 


#Adding a new column to label values to make it easier for plotting
KIRC_FC_P$DifferentiallyExpressed <-"NO Difference"  

#Setting the threshold for p-value and Log fold change 
KIRC_FC_P$DifferentiallyExpressed[KIRC_FC_P$FC > 0.585 & KIRC_FC_P$rawPvalue < 0.05] <- "Upregulated"
KIRC_FC_P$DifferentiallyExpressed[KIRC_FC_P$FC < -0.585 & KIRC_FC_P$rawPvalue < 0.05] <- "Downregulated"

#Using ggplot for plotting and ggrepel for labeling
ggplot(data=KIRC_FC_P, aes(x=FC, y=-log10(rawPvalue), col=DifferentiallyExpressed, label=HUGO )) +
  geom_point()+                      #Creating Scatter plot
  theme_minimal()+                  #Theme with no background annotaions 
  geom_text_repel()+                #Used to avoid overlapping of labels 
  geom_label_repel()+
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-0.585, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")+
  coord_cartesian(xlim =c(-7.5, 7.5), ylim = c(0, 30))+
  ggtitle("KIRC Volcano plot")      #Labeling the Figure 
