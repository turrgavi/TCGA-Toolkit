#Script dump for analysis functions using TCGA gene sets
library("stringr")
library("ggplot2")

#Sample Type Selection function
st_select <- function(sample_type){
  
 
  
  if (sample_type == "Primary" || sample_type == "primary"){
    data <<- subset(data, str_sub(data[,1], -3) == "-01")
    print("Subsetted Primary Solid Tumour Only")
  }
  
  if (sample_type == "Normal" || sample_type == "normal"){
    data <<- subset(data, str_sub(data[,1], -3) == "-11")
  }
  
  if (sample_type == "Metastatic" || sample_type == "metastatic"){
    data <<- subset(data, str_sub(data[,1], -3) == "-06")
  }
  
  if (sample_type == "All" || sample_type == "all"){
    print("Samples have not been subsetted")
  }
  
  if (sample_type == "no Primary" || sample_type == "no primary"){
    print("Primary Tumour samples have been excluded")
    data <<- subset(data, str_sub(data[,1], -3) != "-01")
  }
  
  if (sample_type == "no Metastatic" || sample_type == "no metastatic"){
    print("Metastatic samples have been excluded")
    data <<- subset(data, str_sub(data[,1], -3) != "-06")
  }
  
  if (sample_type == "no Normal" || sample_type == "no normal"){
    print("Normal samples have been excluded")
    data <<- subset(data, str_sub(data[,1], -3) != "-11")
  }
  
}

#GQPlot (Gene Quartile Plot) is used to plot and calculate the P-value of a genes expression based on the Quartile expression of another gene
#Qgene is the gene you wish to subset into upper and lower quartile expression 
#Pgene is the gene you wish to plot seperated by Qgene quartiles
GQPlot <- function(Qgene, Pgene){
  
  Hi_Quart <<- quantile(as.numeric(data[,Qgene]), 0.75)
  Lo_Quart <<- quantile(as.numeric(data[,Qgene]), 0.25)
  
  quartlist <<- vector(length = length(data[,1]))
  x= 1
  for (i in 1:length(data[,1])){
    if (data[i,Qgene] > Hi_Quart){
      quartlist[[x]] <<- "Upper Quartile"
      
      print("UQ")
    }
    
    if (data[i,Qgene] < Lo_Quart){
      quartlist[[x]] <<- "Lower Quartile"
      
      print("LQ")
    }
    x=x+1
    ++i
    print(i)
  }
  
  defined <- as.data.frame(cbind(data[,Pgene], quartlist))
  defined$quartlist <- as.factor(defined$quartlist)
  defined$V1 <- as.numeric(defined$V1)
  defined <- subset(defined, quartlist != FALSE)
  
  ggboxplot(defined, x='quartlist', y='V1', add = "jitter", color = "quartlist", title = Qgene, xlab = paste(Qgene, "Quartile"), ylab = "Normalized Expression Level", palette = c("#66b8fa", "#fa5555")) +
    stat_compare_means(method = "t.test") +
    ggtitle(paste("Expression of", Pgene, "by", Qgene, "Expression Quartile"))+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(color = paste("Quartile"))
  
}

boxplot_gene <- function(gene){
  COX_defined <- as.data.frame(cbind(data[,gene], quartlist))
  COX_defined$quartlist <- as.factor(COX_defined$quartlist)
  COX_defined$V1 <- as.numeric(COX_defined$V1)
  COX_defined <- subset(COX_defined, quartlist != FALSE)
  
  ggboxplot(COX_defined, x='quartlist', y='V1', add="jitter", color = "quartlist")
  
}

#GQPlot but each gene against every other and store p-values
GQPlotAll <- function(){
  cat('FUNCTION NOT WORKING YET')
  for (i in 1:length[,1]){
  Hi_Quart <<- quantile(as.numeric(data[,Qgene]), 0.75)
  Lo_Quart <<- quantile(as.numeric(data[,Qgene]), 0.25)
  
  quartlist <<- vector(length = length(data[,1]))
  x= 1
  for (i in 1:length(data[,1])){
    if (data[i,Qgene] > Hi_Quart){
      quartlist[[x]] <<- "Upper Quartile"
      
      print("UQ")
    }
    
    if (data[i,Qgene] < Lo_Quart){
      quartlist[[x]] <<- "Lower Quartile"
      
      print("LQ")
    }
    x=x+1
    ++i
    print(i)
  }
  
  defined <- as.data.frame(cbind(data[,Pgene], quartlist))
  defined$quartlist <- as.factor(defined$quartlist)
  defined$V1 <- as.numeric(defined$V1)
  defined <- subset(defined, quartlist != FALSE)
  
  ggboxplot(defined, x='quartlist', y='V1', add = "jitter", color = "quartlist", title = Qgene, xlab = paste(Qgene, "Quartile"), ylab = "Normalized Expression Level", palette = c("#66b8fa", "#fa5555")) +
    stat_compare_means(method = "t.test") +
    ggtitle(paste("Expression of", Pgene, "by", Qgene, "Expression Quartile"))+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(color = paste("Quartile"))
  
}

#boxplot_gene <- function(gene){
  #COX_defined <- as.data.frame(cbind(data[,gene], quartlist))
 # COX_defined$quartlist <- as.factor(COX_defined$quartlist)
 # COX_defined$V1 <- as.numeric(COX_defined$V1)
 # COX_defined <- subset(COX_defined, quartlist != FALSE)
  
 # ggboxplot(COX_defined, x='quartlist', y='V1', add="jitter", color = "quartlist")
#}
}

sub_clinical <- function(gene, clin_var){
  
  clin_data <<- read.csv("Clinical_data.csv", header = TRUE)
  
  #clin_data$submitter_id <<- str_sub(clin_data$submitter_id, end = -2)
  
  clin_data <<- subset(clin_data, clin_data[,clin_var] != "")
  
  matched_clin <<- cbind(clin_data$sampleID, clin_data[,clin_var])
  colnames(matched_clin) <- c("sample", "clinical_var")
  
  sub_gene <- data[,gene]
  sub_gene <- as.data.frame(cbind( data[,1], sub_gene))
  colnames(sub_gene) <- c("sample", "gene")
  
  merged <<- merge(sub_gene, matched_clin, by = "sample")
  merged$clinical_var <- as.factor(merged$clinical_var)
  merged$gene <- as.numeric(merged$gene)
  
  variable_count <<- length(unique(merged[,"clinical_var"]))
  if (variable_count >2){
    stat_method <<- "anova"
  }
  
  if (variable_count <= 2){
    stat_method <<- "t.test"
  }
  
  b <- ggboxplot(merged, x='clinical_var', y='gene', 
            add = "jitter", 
            fill = 'clinical_var', 
            title = paste(gene, "expression by clin_var"), 
            xlab = clin_var, 
            ylab = paste("Normalized", 
            gene, "Expression Level")) +
            #palette =c("#139bb3", "#f94343")) +
    
    
    stat_compare_means(method = stat_method, size = 10, label.y = 2.5, label.x = 2 )+
    ggtitle(paste(gene, "Expression by", clin_var))
    
  b +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, size = 36,face="bold"), 
          axis.text=element_text(size=25),
          axis.title=element_text(size=36, margin = margin(r = 20)),
          legend.position = "none")+
    labs(color = paste(clin_var))
}

stats_test <- function(gene){
  pairing <- readline("Is your data paired (P) or unpaired (UP)?: ")
  
  swt_pval <- shapiro.test(as.numeric(data[,gene]))$p.value
  
  if (pairing == "P" && swt_pval >= 0.05){
    method <- "t.test"
    print("Paired T test used")
  }
  
  if (pairing == "UP" && swt_pval >= 0.05){
    method <- "t.test" #FIX THIS
    print("Unpaired T test used")
  }
  
  if (pairing == "P" && swt_pval < 0.05){
    method <- "wilcox"
    print("Wilcox test used")
  }
  
  if (pairing == "UP" && swt_pval < 0.05){
    method <- "mann whitney" #fix this
    print("Mann Whitney test used")
  }
}


#Combine Stage I and II + Stage III and IV
stage_combine_1 <- function(gene){
  clin_data <- read.csv("Clinical_data.csv", header = TRUE)
  
  clin_data$submitter_id <- str_sub(clin_data$submitter_id, end = -2)
  
  run_length <<- length(clin_data[,"clinical_stage"])
  x = 1
  for (x in 1:run_length){
    
    if (clin_data[x,"clinical_stage"] == "Stage I" || clin_data[x,"clinical_stage"] == "Stage II"){
      clin_data[x, "clinical_stage"] <- "(Stage 1 and 2)"
    }
    
    if (clin_data[x,"clinical_stage"] == "Stage III" || clin_data[x,"clinical_stage"] == "Stage IVA" || clin_data[x,"clinical_stage"] == "Stage IVB" || clin_data[x,"clinical_stage"] == "Stage IVC" || clin_data[x,"clinical_stage"] == "Stage IV"){
      clin_data[x, "clinical_stage"] <- "(Stage 3 and 4)"
    }
    
    if (str_sub(clin_data[x,"submitter_id"], -3) == "-11") {
      clin_data[x, "clinical_stage"] <- "Normal Tissue"
      print("normal")
    }
    
  }
  ++x
  
  clin_data <- subset(clin_data, clin_data[,"clinical_stage"] != "")
  
  matched_clin <<- cbind(clin_data$submitter_id, clin_data[,"clinical_stage"])
  colnames(matched_clin) <- c("sample", "clinical_stage")
  
  sub_gene <- data[,gene]
  sub_gene <- as.data.frame(cbind( data[,1], sub_gene))
  colnames(sub_gene) <- c("sample", "gene")
  
  merged <<- merge(sub_gene, matched_clin, by = "sample")
  merged$clinical_stage <- as.factor(merged$clinical_stage)
  merged$gene <- as.numeric(merged$gene)
  
  variable_count <<- length(unique(merged[,"clinical_stage"]))
  if (variable_count >2){
    stat_method <<- "anova"
  }
  
  if (variable_count <= 2){
    stat_method <<- "t.test"
  }
  
  b <- ggboxplot(merged, x='clinical_stage', y='gene', 
                 add = "jitter", 
                 fill = 'clinical_stage', 
                 title = paste(gene, "expression by HPV Status"), 
                 xlab = "Clinical Stage", 
                 ylab = paste("Normalized", 
                              gene, "Expression Level")) +
    #palette =c("#139bb3", "#f94343")) +
    
    
    stat_compare_means(method = stat_method, size = 10, label.y = 3, label.x = 1 )+
    ggtitle(paste(gene, "Expression by Clinical Stage"))
  
  b +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, size = 36,face="bold"), 
          axis.text=element_text(size=25),
          axis.title=element_text(size=36, margin = margin(r = 20)),
          legend.position = "none")+
    labs(color = paste("Clinical Stage"))
  
}

#Combine stage IVs
stage_combine_2 <- function(gene){
  
  clin_data <- read.csv("Clinical_data.csv", header = TRUE)
  
  clin_data$submitter_id <- str_sub(clin_data$submitter_id, end = -2)
  
  run_length <<- length(clin_data[,"clinical_stage"])
  x = 1
  for (x in 1:run_length){
    
    if (clin_data[x,"clinical_stage"] == "Stage IVA" || clin_data[x,"clinical_stage"] == "Stage IVB" || clin_data[x,"clinical_stage"] == "Stage IVC" || clin_data[x,"clinical_stage"] == "Stage IV"){
      clin_data[x, "clinical_stage"] <- "Stage IV"
    }
    
    if (str_sub(clin_data[x,"submitter_id"], -3) == "-11") {
      clin_data[x, "clinical_stage"] <- "Normal Tissue"
      print("normal")
    }
  }
  ++x
  
  clin_data <- subset(clin_data, clin_data[,"clinical_stage"] != "")
  
  matched_clin <<- cbind(clin_data$submitter_id, clin_data[,"clinical_stage"])
  colnames(matched_clin) <- c("sample", "clinical_stage")
  
  sub_gene <- data[,gene]
  sub_gene <- as.data.frame(cbind( data[,1], sub_gene))
  colnames(sub_gene) <- c("sample", "gene")
  
  merged <<- merge(sub_gene, matched_clin, by = "sample")
  merged$clinical_stage <- as.factor(merged$clinical_stage)
  merged$gene <- as.numeric(merged$gene)

  variable_count <<- length(unique(merged[,"clinical_stage"]))
  if (variable_count >2){
    stat_method <<- "anova"
  }
  
  if (variable_count <= 2){
    stat_method <<- "t.test"
  }
  
  b <- ggboxplot(merged, x='clinical_stage', y='gene', 
                 add = "jitter", 
                 fill = 'clinical_stage', 
                 title = paste(gene, "expression by HPV Status"), 
                 xlab = "Clinical Stage", 
                 ylab = paste("Normalized", 
                              gene, "Expression Level")) +
    #palette =c("#139bb3", "#f94343")) +
    
    
    stat_compare_means(method = stat_method, size = 10, label.y = 3, label.x = 1)+
    ggtitle(paste(gene, "Expression by Clinical Stage"))
  
  b +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, size = 36,face="bold"), 
          axis.text=element_text(size=25),
          axis.title=element_text(size=36, margin = margin(r = 20)),
          legend.position = "none")+
    labs(color = paste("Clinical Stage"))
 

}

gene_corr <- function(cor_gene1, cor_gene2){
  corr_data <- as.data.frame(cbind(as.numeric(data[,cor_gene1]), as.numeric(data[,cor_gene2])))
  
  colnames(corr_data) <- c(cor_gene1, cor_gene2)
  rownames(corr_data) <- data[,1]
  
  
  corr_value <- cor.test(corr_data[,cor_gene1], corr_data[,cor_gene2])
  ggscatter(corr_data, x = cor_gene1, y = cor_gene2, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = cor_gene1, ylab = cor_gene2,
            ggtheme = theme_minimal_grid())
}

