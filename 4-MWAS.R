# HEADER -----------------------------------------------------------------------
#
# TITLE:   HMOxmiRNA.R
#
# PURPOSE: Estimate associations between breastmilk miRNA and HMOs
#
# INPUT:   miRNA_cpm.csv
#          meta_clean_bl.csv
#
# OUTPUT:  
#
# DATE:    March 24, 2022
#
#
#          
# SET UP -----------------------------------------------------------------------
#clear the workspace
rm(list = ls())

#turn off scientific notation
options(scipen = 100)

#load libraries
library(plyr); library(tidyverse); library(purrr); library(readxl);
library(stringr); library(lme4); library(lmerTest); library(corrplot); 
library(lubridate); library(MASS); library(pscl); library(ggbiplot)
library(stargazer); library(ggrepel); library(openxlsx)

#set the input folder
data_in <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/rda"

#set the output folder
figs_out <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/figs"

#read in the clean meta data
meta <- read.csv(file = paste0(data_in, "meta_clean_EY.csv"))
meta <- read.csv("input/meta_clean_EY.csv")

#miRNA_cpm
miRNA_cpm <- read.csv(file = paste0(data_in, "miRNA_counts_EY.csv"))
miRNA_cpm <- read.csv("input/miRNA_counts_EY.csv")

# FORMAT THE DATA --------------------------------------------------------------

#drop the extra first column
miRNA_cpm <- miRNA_cpm[,-1]

#make a new data frame that excludes the QC variables
#add an X to row names so R plays nice
row.names(miRNA_cpm) <- paste0("X", miRNA_cpm$dyad_id)
miRNA_cpm_noQC <- miRNA_cpm[,1:210]

#drop the X column 
meta <- meta[,-1]
#update the format of dyad_id
meta$dyad_id <- as.numeric(substr(meta$merge_id_dyad,4,7))

#combine the miRNA data with the meta data
meta_miRNA <- left_join(miRNA_cpm, meta, by = "dyad_id")

# make a list of summary vars and HMOs concentrations to loop over
HMOs <- c("Secretor", "Diversity", "Sia", "Fuc", "x2FL_nmol_ml","x3FL_nmol_ml",
          "LNnT_nmol_ml",
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

HMOs <- unlist(HMOs)

# SECRETOR REGRESSION ----------------------------------------------------------
# make a data frame to store the results
table(meta_miRNA$Secretor)
# No Yes 
# 12  96
sec_result <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(sec_result) <- c("est", "CI_lower", "CI_upper", "p", "miRNA")

pdf("figs/plots_secretor.pdf")

for(i in 1:210){ # loop over each miRNA
  tryCatch(
    {
      fit <- lm(meta_miRNA[,i] ~ meta_miRNA$Secretor + prop_rRNA + 
                  Vol_Supernatant + Date_Evs, data = meta_miRNA)
      #save the results to the ith row in the output table
      #estimate
      sec_result[i, 1] <- round(summary(fit)$coefficients[2,1], digits = 2) 
      #lower confidence interval
      sec_result[i, 2] <- round(confint(fit)[2,1], digits = 2)
      #upper confidence interval
      sec_result[i, 3] <- round(confint(fit)[2,2], digits = 2) 
      #p-value
      sec_result[i, 4] <- round(summary(fit)$coefficients[2,4], digits = 10)
      #name of miRNA
      sec_result[i, 5] <- colnames(meta_miRNA)[i]
      
      if(!is.na(sec_result[i, 4]) & sec_result[i,4] < 0.05){
        plot(meta_miRNA[,i] ~ meta_miRNA$Secretor,
             xlab = "Secretor",
             ylab = colnames(meta_miRNA[i]),
             cex.lab = 0.5)
      }
      
      {print(i)}
    }, error=function(e){})
}

dev.off()
warnings()

sum(sec_result$p < 0.05)
# 15 miRNAs were associated with secretor status
# now only 3 

print(sec_result$miRNA[which(sec_result$p < 0.05)])
#"hsa.miR.152.3p"  "hsa.miR.361.3p"  "hsa.miR.1287.5p"

#adjust for multiple testing
sec_result$p_adj <- p.adjust(sec_result$p, method = "BH")
sum(sec_result$p_adj < 0.05)
#NA

# HMO REGRESSION ---------------------------------------------

# make an Excel workbook to store the results
#wb <- createWorkbook(type = "xlsx")
wb <- createWorkbook()
# Save the workbook as an XLSX file
#saveWorkbook(wb, "your_filename.xlsx")

# run a loop to perform MWAS on each HMO
for(thisHMO in HMOs){ # loop over each HMO
  
  # make a data frame to store results
  thisResult <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(thisResult) <- c("est", "CI_lower", "CI_upper", "p", "miRNA")
  
  # make a pdf for plots
  pdf("hmo_mwas_plots", thisHMO, ".pdf")
  
  for(i in 1:210){ # loop over each miRNA
    tryCatch(
      {
        fit <- lm(meta_miRNA[,thisHMO] ~ meta_miRNA[,i] + prop_rRNA + 
                    Vol_Supernatant + Date_Evs + Secretor, data = meta_miRNA)
        #save the results to the ith row in the output table
        #estimate
        thisResult[i, 1] <- round(summary(fit)$coefficients[2,1], digits = 2) 
        #lower confidence interval
        thisResult[i, 2] <- round(confint(fit)[2,1], digits = 2)
        #upper confidence interval
        thisResult[i, 3] <- round(confint(fit)[2,2], digits = 2) 
        #p-value
        thisResult[i, 4] <- round(summary(fit)$coefficients[2,4], digits = 10)
        #name of miRNA
        thisResult[i, 5] <- colnames(meta_miRNA)[i]
        
        if(!is.na(thisResult[i, 4]) & thisResult[i,4] < 0.05){
          plot(meta_miRNA[,thisHMO] ~ meta_miRNA[,i],
               xlab = colnames(meta_miRNA[i]),
               ylab = thisHMO,
               cex.lab = 0.5)
        }
        
        {print(i)}
      }, error=function(e){})
  }
  
  dev.off()
  warnings()
  
  # calculate the corrected p value
  thisResult$pfdr <- round(p.adjust(thisResult$p, method = "BH", 
                                    n = length(thisResult$p)*20), digits = 10) 
  # multiply length of thisResult by 20 to account for the 20 HMOs
  #sort by pfdr value
  thisResult <- thisResult[order(thisResult$pfdr),]
  
  
  #save the output to a new sheet in wb
  # sheet <- createSheet(wb, sheetName = thisHMO)
  sheet <- addWorksheet(wb, sheetName = thisHMO)
  #addDataFrame(thisResult, sheet)
  writeData(wb, sheet, thisResult, startRow = 1, 
            startCol = 1, rowNames = FALSE)
  #write.xlsx(thisResult, sheet = sheet, startRow = 1, 
  #          startCol = 1, rowNames = FALSE)
  
  # if the HMO is one that we want to plot, save the results into a data frame
  if(thisHMO == "X2FL_ug_ml"){
    X2FL_out <- data.frame(thisResult)
  }
  if(thisHMO == "X3FL_ug_ml"){
    X3FL_out <- data.frame(thisResult)
  }
  if(thisHMO == "X3SL_ug_ml"){
    X3SL_out <- data.frame(thisResult)
  }
  if(thisHMO == "DSLNT_nmol_ml") {
    DSLNT_out <- data.frame(thisResult)
  }

}

# save the results to Excel (commented out to avoid unintentional overwriting)
# saveWorkbook(wb, paste0(figs_out,"HMO_MWAS_results.xlsx"))
# HMO_MWAS_results.xlsx ####

# Volcano plots ----------------------------------------------------------------

# make volcano plots for 2'fucosyllactose, 3'-fucosyllactose, 3'sialyllactose 
# and DSLNT
# drop the lines that didn't converge
X2FL_out <- X2FL_out[complete.cases(X2FL_out),]

# add variables to make the color 
X2FL_out$Direction <- NA
X2FL_out$Direction[(X2FL_out$est > 0 & X2FL_out$pfdr < 0.05)] <- "up"
X2FL_out$Direction[(X2FL_out$est < 0 & X2FL_out$pfdr < 0.05)] <- "down"

X2FL_out$label <- NA
X2FL_out$label[(X2FL_out$pfdr < 0.05)] <- substr(X2FL_out$miRNA,5,20)[X2FL_out$pfdr < 0.05]

# plot
plot_X2FL_out <- ggplot(data = X2FL_out, aes(x = est, y = -log10(pfdr), 
                                               col = Direction,
                                               label = label)) + 
  ylab("-log10(pfdr)") + 
  xlab("Fold Change") + 
  geom_point() +
  geom_text_repel(data = X2FL_out, mapping = aes(label = label)) +
  theme_minimal() +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("blue", "red", "black")) +
  theme(legend.position = "none")

plot_X2FL_out

# 3' fucosyllactose

X3FL_out <- X3FL_out[complete.cases(X3FL_out),]

# add variables to make the color 
X3FL_out$Direction <- NA
X3FL_out$Direction[(X3FL_out$est > 0 & X3FL_out$pfdr < 0.05)] <- "up"
X3FL_out$Direction[(X3FL_out$est < 0 & X3FL_out$pfdr < 0.05)] <- "down"

X3FL_out$label <- NA
X3FL_out$label[(X3FL_out$pfdr < 0.05)] <- substr(X3FL_out$miRNA,5,20)[X3FL_out$pfdr < 0.05]

# plot
plot_X3FL_out <- ggplot(data = X3FL_out, aes(x = est, y = -log10(pfdr), 
                                             col = Direction,
                                             label = label)) + 
  ylab("-log10(pfdr)") + 
  xlab("Fold Change") + 
  geom_point() +
  geom_text_repel(data = X3FL_out, mapping = aes(label = label)) +
  theme_minimal() +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("blue", "red", "black")) +
  theme(legend.position = "none")

plot_X3FL_out

# X3SL

X3SL_out <- X3SL_out[complete.cases(X3SL_out),]

# add variables to make the color 
X3SL_out$Direction <- NA
X3SL_out$Direction[(X3SL_out$est > 0 & X3SL_out$pfdr < 0.05)] <- "up"
X3SL_out$Direction[(X3SL_out$est < 0 & X3SL_out$pfdr < 0.05)] <- "down"

X3SL_out$label <- NA
X3SL_out$label[(X3SL_out$pfdr < 0.05)] <- substr(X3SL_out$miRNA,5,20)[X3SL_out$pfdr < 0.05]

# plot
plot_X3SL_out <- ggplot(data = X3SL_out, aes(x = est, y = -log10(pfdr), 
                                             col = Direction,
                                             label = label)) + 
  ylab("-log10(pfdr)") + 
  xlab("Fold Change") + 
  geom_point() +
  geom_text_repel(data = X3SL_out, mapping = aes(label = label)) +
  theme_minimal() +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("blue", "red", "black")) +
  theme(legend.position = "none")

plot_X3SL_out

# X3SL

DSLNT_out <- DSLNT_out[complete.cases(DSLNT_out),]

# add variables to make the color 
DSLNT_out$Direction <- NA
DSLNT_out$Direction[(DSLNT_out$est > 0 & DSLNT_out$pfdr < 0.05)] <- "up"
DSLNT_out$Direction[(DSLNT_out$est < 0 & DSLNT_out$pfdr < 0.05)] <- "down"

DSLNT_out$label <- NA
DSLNT_out$label[(DSLNT_out$pfdr < 0.05)] <- substr(DSLNT_out$miRNA,5,20)[DSLNT_out$pfdr < 0.05]

# plot
plot_DSLNT_out <- ggplot(data = DSLNT_out, aes(x = est, y = -log10(pfdr), 
                                             col = Direction,
                                             label = label)) + 
  ylab("-log10(pfdr)") + 
  xlab("Fold Change") + 
  geom_point() +
  geom_text_repel(data = DSLNT_out, mapping = aes(label = label)) +
  theme_minimal() +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("blue", "red", "black")) +
  theme(legend.position = "none")

plot_DSLNT_out

# REMOVE NON-SECRETORS ---------------------------------------------------------

# make an Excel workbook to store the results
wb_sec <- createWorkbook(type = "xlsx")

meta_sec <- meta_miRNA[which(meta_miRNA$Secretor == "Yes"),]
# run a loop to perform MWAS on each HMO
for(thisHMO in HMOs){ # loop over each HMO
  
  # make a data frame to store results
  thisResult <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(thisResult) <- c("est", "CI_lower", "CI_upper", "p", "miRNA")
  
  # make a pdf for plots
  pdf(paste0(figs_out, "plots_", thisHMO, ".pdf"))
  
  for(i in 1:210){ # loop over each miRNA
    tryCatch(
      {
        fit <- lm(meta_sec[,thisHMO] ~ meta_sec[,i] + prop_rRNA + 
                    Vol_Supernatant + Date_Evs, data = meta_sec)
        #save the results to the ith row in the output table
        #estimate
        thisResult[i, 1] <- round(summary(fit)$coefficients[2,1], digits = 2) 
        #lower confidence interval
        thisResult[i, 2] <- round(confint(fit)[2,1], digits = 2)
        #upper confidence interval
        thisResult[i, 3] <- round(confint(fit)[2,2], digits = 2) 
        #p-value
        thisResult[i, 4] <- round(summary(fit)$coefficients[2,4], digits = 10)
        #name of miRNA
        thisResult[i, 5] <- colnames(meta_miRNA)[i]
        
        if(!is.na(thisResult[i, 4]) & thisResult[i,4] < 0.05){
          plot(meta_sec[,thisHMO] ~ meta_sec[,i],
               xlab = colnames(meta_sec[i]),
               ylab = thisHMO,
               cex.lab = 0.5)
        }
        
        {print(i)}
      }, error=function(e){})
  }
  
  dev.off()
  warnings()
  
  # calculate the corrected p value
  thisResult$pfdr <- round(p.adjust(thisResult$p, method = "BH", 
                                    n = length(thisResult$p)*19), digits = 10) 
  # multiply length of thisResult by 20 to account for the 19 HMOs
  #sort by pfdr value
  thisResult <- thisResult[order(thisResult$pfdr),]
  
  
  #save the output to a new sheet in wb
  sheet <- createSheet(wb_sec, sheetName = thisHMO)
  addDataFrame(thisResult, sheet)
  
  # if the HMO is one that we want to plot, save the results into a data frame
  if(thisHMO == "X2FL_ug_ml"){
    X2FL_out <- data.frame(thisResult)
  }
  if(thisHMO == "X3FL_ug_ml"){
    X3FL_out <- data.frame(thisResult)
  }
  if(thisHMO == "X3SL_ug_ml"){
    X3SL_out <- data.frame(thisResult)
  }
  if(thisHMO == "DSLNT_nmol_ml") {
    DSLNT_out <- data.frame(thisResult)
  }
  
}

# save the results to Excel (commented out to avoid unintentional overwriting)
# saveWorkbook(wb_sec, paste0(figs_out,"HMO_MWAS_results_secretors.xlsx"))
# HMO_MWAS_results_secretors.xlsx ####

# REGRESSION, ADJUSTED ---------------------------------------------------------

# run a loop to perform MWAS on each HMO
for(thisHMO in HMOs){ # loop over each HMO
  
  # make a data frame to store results
  thisResult <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(thisResult) <- c("est", "CI_lower", "CI_upper", "p", "miRNA")
  
  # make a pdf for plots
  pdf(paste0(figs_out, "plots_", thisHMO, ".pdf"))
  
  for(i in 1:210){ # loop over each miRNA
    tryCatch(
      {
        fit <- lm(meta_miRNA[,thisHMO] ~ meta_miRNA[,i] + prop_rRNA + 
                    Vol_Supernatant + Date_Evs + Secretor + 
                    breastfeedings_continuous + breast_milk_time_hrs, 
                  data = meta_miRNA)
        #save the results to the ith row in the output table
        #estimate
        thisResult[i, 1] <- round(summary(fit)$coefficients[2,1], digits = 2) 
        #lower confidence interval
        thisResult[i, 2] <- round(confint(fit)[2,1], digits = 2)
        #upper confidence interval
        thisResult[i, 3] <- round(confint(fit)[2,2], digits = 2) 
        #p-value
        thisResult[i, 4] <- round(summary(fit)$coefficients[2,4], digits = 10)
        #name of miRNA
        thisResult[i, 5] <- colnames(meta_miRNA)[i]
        
        if(!is.na(thisResult[i, 4]) & thisResult[i,4] < 0.05){
          plot(meta_miRNA[,thisHMO] ~ meta_miRNA[,i],
               xlab = colnames(meta_miRNA[i]),
               ylab = thisHMO,
               cex.lab = 0.5)
        }
        
        {print(i)}
      }, error=function(e){})
  }
  
  dev.off()
  warnings()
  
  # calculate the corrected p value
  thisResult$pfdr <- round(p.adjust(thisResult$p, method = "BH", 
                                    n = length(thisResult$p)*20), digits = 10) 
  # multiply length of thisResult by 20 to account for the 20 HMOs
  #sort by pfdr value
  thisResult <- thisResult[order(thisResult$pfdr),]
  
  
  #save the output to a new sheet in wb
  sheet <- createSheet(wb, sheetName = thisHMO)
  addDataFrame(thisResult, sheet)
  
  # if the HMO is one that we want to plot, save the results into a data frame
  if(thisHMO == "X3FL..nmol.mL."){
    X3FL_out <- data.frame(thisResult)
  }
  if(thisHMO == "X3.SL..nmol.mL."){
    X3SL_out <- data.frame(thisResult)
  }
}
