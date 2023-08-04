#Import all required packages
library("stringr")
library("ggplot2")
library(cowplot)
library(ggpubr)
source("C:/Users/gavin/OneDrive/University of Queensland/PhD/RScripts/TCGA_Functions.R")

#HNC
#setwd("C:/Users/gavin/OneDrive/University of Queensland/Research/Dr Janin Chandra/TCGA Head and Neck")

#CESC
setwd("C:/Users/gavin/OneDrive/University of Queensland/Research/Dr Janin Chandra/TCGA Cervical Cancer/R")

#LUSC
#setwd("C:/Users/gavin/OneDrive/University of Queensland/Research/Dr Janin Chandra/TCGA Lung Squamous Cell")

#Open file prompt to select raw expression data file
#Organize Matrix
print("Select expression data file as CSV")
raw_exp <- file.choose()
data <- as.matrix(t(read.csv(raw_exp, header = FALSE)))

colnames(data) <- data[1,]
data <- data[-1,]

#rownames(data) <- data[,1]
#data <- data[,-1]

#Subset only primary tumours
#data <- subset(data, str_sub(data[,1], -3) == "-01")

#subset sample type
st_select("primary")


#Plot gene according to another genes expression quartiles
#GQPlot("PTGS2", "ITGB6")

sub_clinical("AURKA", "sample_type")

stage_combine_1("MYC")

stage_combine_2("MYC")

gene_corr("TGFB1", "MYC")
