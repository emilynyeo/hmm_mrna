# HEADER -----------------------------------------------------------------------
#
# TITLE:   2-PCA.R
#
# PURPOSE: do a PCA analysis
#
# DATE:    March 14, 2022
#
# SET UP -----------------------------------------------------------------------
#clear the workspace
rm(list = ls())

#turn off scientific notation
options(scipen = 100)

#load libraries
install.packages('devtools')
library(devtools)
install_github('vqv/ggbiplot')
pacman::p_load(knitr, tidyverse, magrittr, lme4, lmerTest, GGally, corrplot, 
               Hmisc, kableExtra, dplyr, plyr, janitor, lubridate, survminer, 
               ggplot2, here, readr, tableone, officer, flextable,finalfit,
               purrr, stringr, lme4, corrplot, pscl, stargazer, MASS, lmerTest,
               readxl, factoextra, ggbiplot)
#Error in library(qiime2R) : there is no package called ‘qiime2R’

#set the input folder
#data_in <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/rda"
data_in <- "/input/"

#set the output folder
#figs_out <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/figs"
figs_out <- "/figs/"

#miRNA_cpm####
#miRNA_cpm <- read.csv(file = paste0(data_in, "miRNA_count.csv"))
miRNA_cpm_old <- read.csv("/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/rdamiRNA_counts.csv")
miRNA_cpm <- read.csv("input/miRNA_counts_EY.csv")

#drop the extra first column
miRNA_cpm <- miRNA_cpm[,-1]

#make a new data frame that excludes the QC variables
#add an X to row names so R plays nice
row.names(miRNA_cpm) <- paste0("X", miRNA_cpm$dyad_id)
miRNA_cpm_noQC <- miRNA_cpm[,1:210]

# calculate principle components
# miRNA.pca <- prcomp(miRNA_cpm_noQC, center = T, scale = T) # non numeric error
miRNA.pca <- prcomp(miRNA_cpm_noQC[,2:210], center = T, scale = T)

# check whether we can use eigenvalue > 1 as a metric
I((summary(miRNA.pca)$sdev)^2)
# eigenvalue = standard deviation^2. Too many to use those > 1.
plot(miRNA.pca)

temp <- data.frame(t(summary(miRNA.pca)$importance))
temp$PC <- row.names(temp)
temp <- temp[1:10,]
temp$PC <- factor(temp$PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                      "PC7", "PC8", "PC9", "PC10"))

plot_PC <- ggplot(data = temp, aes(x = PC, y = Proportion.of.Variance)) + 
  geom_bar(stat = "identity") +
  ylab("Proportion of Variance") +
  xlab("")
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank()) +
  theme(legend.position = "none")

plot_PC

png("figs/Supplemental_Figure_1_EY.png", 
    units = "in",
    width = 10,
    height = 8,
    pointsize = 12,
    res = 1200)

plot_PC

dev.off()
# Supplemental_Figure_1.png####

# get the loading factors for pca
fviz_pca_var(miRNA.pca, col.var = "contrib",
             select.var = list(contrib = 5),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T)

res.var <- get_pca_var(miRNA.pca)

#save the PCA results in a data frame
PCA <- data.frame(miRNA.pca$x)

#read in the clean meta data
meta <- read.csv("input/meta_clean_EY.csv", stringsAsFactors = TRUE)

# force 'On Time' gestational age to be the reference
summary(meta$gestational_age_cat)
meta$gestational_age_cat <- relevel(meta$gestational_age_cat, 
                                    ref = "On Time")
meta$pred_bf <- relevel(meta$pred_bf, ref = "Yes")

#drop the X column 
meta <- meta[,-1]
#add an X to dyad_id so it merges with PCA
meta$dyad_id <- paste0("X", as.numeric(substr(meta$merge_id_dyad,4,7)))

#only keep individuals in meta who have miRNA data
meta <- meta[meta$dyad_id %in% row.names(PCA),]

PCA$dyad_id <- row.names(PCA)

PCA <- dplyr::select(PCA, c("dyad_id", "PC1":"PC5"))

meta_miRNA <- left_join(meta, PCA, by = "dyad_id")

miRNA_cpm$dyad_id <- paste0("X", miRNA_cpm$dyad_id)

#meta_miRNA <- left_join(meta_miRNA, miRNA_cpm[211:217], by = "dyad_id")
#editting the line above:
miRNA_cpm_cols <- dplyr::select(miRNA_cpm, c("dyad_id","LibSize", 
                                    "LibSizeNormalized", "file_name","LL", 
                                    "Address","Barcode",
                                    "prop_rRNA","Vol_Supernatant","Date_Evs"))
meta_miRNA <- left_join(meta_miRNA, miRNA_cpm_cols, by = "dyad_id")

colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC1"] <- "PC1"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC2"] <- "PC2"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC3"] <- "PC3"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC4"] <- "PC4"
colnames(meta_miRNA)[colnames(meta_miRNA) == "PCA$PC5"] <- "PC5"

# REGRESS PC1 AND PC2 AGAINST HMO SUMMARY VARS ---------------------------------

# compare PCs across secretor status
fit <- lm(meta_miRNA$PC1 ~ meta_miRNA$Secretor + prop_rRNA + Vol_Supernatant + 
            Date_Evs, data = meta_miRNA)
summary(fit)
fit2 <- lm (meta_miRNA$PC2 ~ meta_miRNA$Secretor + prop_rRNA + Vol_Supernatant + 
            Date_Evs, data = meta_miRNA)
summary(fit2)
fit3 <- lm(meta_miRNA$PC3 ~ meta_miRNA$Secretor + prop_rRNA + Vol_Supernatant + 
            Date_Evs, data = meta_miRNA)
summary(fit3)
fit4 <- lm(meta_miRNA$PC4 ~ meta_miRNA$Secretor + prop_rRNA + Vol_Supernatant + 
            Date_Evs, data = meta_miRNA)
summary(fit4)
fit5 <- lm(meta_miRNA$PC5 ~ meta_miRNA$Secretor + prop_rRNA + Vol_Supernatant + 
            Date_Evs, data = meta_miRNA)
summary(fit5)

var_list <- c("Diversity", "Fuc", "Sia")

var_list <- unlist(var_list)

result_unadj <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(result_unadj) <- c("Var", "PC1_Est", "PC1_CI_lower", "PC1_CI_upper", 
                      "PC1_p", "PC2_Est", "PC2_CI_lower", "PC2_CI_upper",
                      "PC2_p")

for(thisVar in var_list){
  
    fit <- lm(meta_miRNA[,thisVar] ~ PC1 + PC2 + prop_rRNA + 
                Vol_Supernatant + Date_Evs + Secretor, data = meta_miRNA)
    #save the results to the next row in the output table
    i <- nrow(result_unadj) + 1
    #variable
    result_unadj[i, 1] <- thisVar
    #estimate PC1
    result_unadj[i, 2] <- round(summary(fit)$coefficients[2,1], digits = 4) 
    #lower confidence interval PC1
    result_unadj[i, 3] <- round(confint(fit)[2,1], digits = 4)
    #upper confidence interval PC1
    result_unadj[i, 4] <- round(confint(fit)[2,2], digits = 4) 
    #p-value PC1
    result_unadj[i, 5] <- round(summary(fit)$coefficients[2,4], digits = 5)
    #estimate PC2
    result_unadj[i, 6] <- round(summary(fit)$coefficients[3,1], digits = 4) 
    #lower confidence interval PC1
    result_unadj[i, 7] <- round(confint(fit)[3,1], digits = 4)
    #upper confidence interval PC1
    result_unadj[i, 8] <- round(confint(fit)[3,2], digits = 4) 
    #p-value PC1
    result_unadj[i, 9] <- round(summary(fit)$coefficients[3,4], digits = 5)
}

result_adj <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(result_adj) <- c("Var", "PC1_Est", "PC1_CI_lower", "PC1_CI_upper", 
                            "PC1_p", "PC2_Est", "PC2_CI_lower", "PC2_CI_upper",
                            "PC2_p")

for(thisVar in var_list){
  
  fit <- lm(meta_miRNA[,thisVar] ~ PC1 + PC2 + prop_rRNA + 
              Vol_Supernatant + Date_Evs + age_in_days + Secretor +
              breast_milk_time_hrs + pred_bf, data = meta_miRNA)
  #save the results to the next row in the output table
  i <- nrow(result_adj) + 1
  #variable
  result_adj[i, 1] <- thisVar
  #estimate PC1
  result_adj[i, 2] <- round(summary(fit)$coefficients[2,1], digits = 4) 
  #lower confidence interval PC1
  result_adj[i, 3] <- round(confint(fit)[2,1], digits = 4)
  #upper confidence interval PC1
  result_adj[i, 4] <- round(confint(fit)[2,2], digits = 4) 
  #p-value PC1
  result_adj[i, 5] <- round(summary(fit)$coefficients[2,4], digits = 5)
  #estimate PC2
  result_adj[i, 6] <- round(summary(fit)$coefficients[3,1], digits = 4) 
  #lower confidence interval PC1
  result_adj[i, 7] <- round(confint(fit)[3,1], digits = 4)
  #upper confidence interval PC1
  result_adj[i, 8] <- round(confint(fit)[3,2], digits = 4) 
  #p-value PC1
  result_adj[i, 9] <- round(summary(fit)$coefficients[3,4], digits = 5)
}

block_table(result_adj)
block_table(result_unadj)

# regress indiidual HMO concentrations against PC1 and PC2
var_list2 <- c("x2FL_nmol_ml","x3FL_nmol_ml","LNnT_nmol_ml",
               "x3SL_nmol_ml","DFLac_nmol_ml",
               "x6SL_nmol_ml","LNT_nmol_ml","LNFP_I_nmol_ml",
               "LNFP_II_nmol_ml","LNFP_III_nmol_ml",
               "LSTb_nmol_ml","LSTc_nmol_ml","DFLNT_nmol_ml","LNH_nmol_ml",
               "DSLNT_nmol_ml","FLNH_nmol_ml","DFLNH_nmol_ml","FDSLNH_nmol_ml",
               "DSLNH_nmol_ml","SUM_nmol_ml",
               "Sia_nmol_ml","Fuc_nmol_ml","x2FL_ug_ml","x3FL_ug_ml","LNnT_ug_ml",
               "x3SL_ug_ml", "DFLac_ug_ml", "x6SL_ug_ml", "LNT_ug_ml", "LNFP_I_ug_ml",
               "LNFP_II_percent","LNFP_III_percent","LSTb_percent","LSTc_percent",
               "DFLNT_percent","LNH_percent", "DSLNT_percent", "FLNH_percent", 
               "DFLNH_percent", "FDSLNH_percent", "DSLNH_percent")

var_list2 <- unlist(var_list)

result_HMO <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(result_HMO) <- c("Var", "PC1_Est", "PC1_CI_lower", "PC1_CI_upper", 
                            "PC1_p", "PC2_Est", "PC2_CI_lower", "PC2_CI_upper",
                            "PC2_p")

for(thisVar in var_list2){
  
  fit <- lm(meta_miRNA[,thisVar] ~ PC1 + PC2 + prop_rRNA + 
              Vol_Supernatant + Date_Evs + Secretor + breast_milk_time_hrs +
              age_in_days + pred_bf, data = meta_miRNA)
  #save the results to the next row in the output table
  i <- nrow(result_HMO) + 1
  #variable
  result_HMO[i, 1] <- thisVar
  #estimate PC1
  result_HMO[i, 2] <- round(summary(fit)$coefficients[2,1], digits = 4) 
  #lower confidence interval PC1
  result_HMO[i, 3] <- round(confint(fit)[2,1], digits = 4)
  #upper confidence interval PC1
  result_HMO[i, 4] <- round(confint(fit)[2,2], digits = 4) 
  #p-value PC1
  result_HMO[i, 5] <- round(summary(fit)$coefficients[2,4], digits = 5)
  #estimate PC2
  result_HMO[i, 6] <- round(summary(fit)$coefficients[3,1], digits = 4) 
  #lower confidence interval PC2
  result_HMO[i, 7] <- round(confint(fit)[3,1], digits = 4)
  #upper confidence interval PC2
  result_HMO[i, 8] <- round(confint(fit)[3,2], digits = 4) 
  #p-value PC2
  result_HMO[i, 9] <- round(summary(fit)$coefficients[3,4], digits = 5)
}

# adjust for multiple testing
result_HMO$PC1_pfdr <- p.adjust(result_HMO$PC1_p, method = "BH")
result_HMO$PC2_pfdr <- p.adjust(result_HMO$PC2_p, method = "BH")

# visualizations for Secretor = Yes vs. No
ggbiplot(miRNA.pca, group = meta_miRNA$Secretor, ellipse = T, 
         var.axes = F) +
  theme_minimal() +
  labs(color = "Breast\nFeeding\nCategory") +
  ylim(-4,2.5) + xlim(-3,3) +
  ylab("PC2 (15.6% explained var.)") + xlab("PC1 (18.6% explained var.)")

# visualizations for breast feedings per day vs formula feedings per day
# fix the order before plotting
meta_miRNA$breastfeedingcat <- factor(meta$breastfeedingcat, 
                                      levels = c("High", "Med", "Low", "None"))
meta_miRNA$formulacat <- factor(meta$formulacat, 
                                levels = c("High", "Med", "Low", "None"))

png(filename = paste0(figs_out, "Figure_1_A_EY.png"), 
    units = "in",
    width = 5,
    height = 4,
    pointsize = 12,
    res = 1200)
ggbiplot(miRNA.pca, group = meta_miRNA$breastfeedingcat, ellipse = T, 
         var.axes = F) +
  theme_minimal() +
  labs(color = "Breast\nFeeding\nCategory") +
  ylim(-4,2.5) + xlim(-3,3) +
  ylab("PC2 (15.6% explained var.)") + xlab("PC1 (18.6% explained var.)")
dev.off()

png(filename = paste0(figs_out, "Figure_1_B_EY.png"), 
    units = "in",
    width = 5,
    height = 4,
    pointsize = 16,
    res = 1200)
ggbiplot(miRNA.pca, group = meta_miRNA$formulacat, ellipse = T, 
         var.axes = F) +
  ylim(-4,2.5) + xlim(-3,3) +
  ylab("PC2 (15.6% explained var.)") + xlab("PC1 (18.6% explained var.)") +
  theme_minimal() +
  labs(color = "Formula\nCategory")
dev.off()
