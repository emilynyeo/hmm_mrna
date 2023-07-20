# HEADER -----------------------------------------------------------------------
#
# TITLE:   1-data-cleaning.R
#
# PURPOSE: clean the meta-data, miRNA_cpm, and raw miRNA counts
#
# DATE:    March 22, 2022
#
# SET UP -----------------------------------------------------------------------

#clear the workspace
rm(list = ls())

#turn off scientific notation
options(scipen = 100)

#load libraries
pacman::p_load(knitr, tidyverse, magrittr, lme4, lmerTest, GGally, corrplot, 
               Hmisc, kableExtra, dplyr, plyr, janitor, lubridate, survminer, 
               ggplot2, here, readr, tableone, officer, flextable,finalfit,
               purrr, stringr, lme4, corrplot, pscl, stargazer, MASS, lmerTest,
               readxl)
#install.packages("rJava", repos="https://rforge.net")
#library(xlsx) #this isn't working

#set up output folders
#rda_out for cleaned data sets
#rda_out <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/rda"
rda_out <- "/ouput"

#figs_out for any figures or tables
#figs_out <- "/Volumes/IPHY/ADORLab/__Users/emye7956/HMO-miRNA/1-data-cleaning/figs"
figs_out <- "/figs"

# read in the miRNA data
#counts per million for use when miRNA is the predictor
miRNA_cpm <- read.csv("/Volumes/IPHY/ADORLab/Lab\ Projects/Mothers\ Milk/Papers/miRNA\ Methods/Input/210914_Alderette_Goran_exceRpt_70perc_expr_miRNAs_TMMnormalized_CPM_MOM155removed.csv")
#QC file contains important control variables
miRNA_QC <- read.csv("/Volumes/IPHY/ADORLab/Lab\ Projects/Mothers\ Milk/Papers/miRNA\ Methods/Input/Goran_Alderette_Pilot_9-2021_Batch_Info_Combined.csv")
#raw counts for use when miRNA is the outcome
miRNA_counts <- read.csv("/Volumes/IPHY/ADORLab/Lab\ Projects/Mothers\ Milk/Papers/miRNA\ Methods/Input/210914_Alderette_Goran_exceRpt_70perc_expr_miRNAs_TMMnormalized_RawCounts_MOM155removed.csv")

# read in the unfiltered miRNA data frame
miRNA_all <- read.csv("/Volumes/IPHY/ADORLab/Lab\ Projects/Mothers\ Milk/Breast\ Milk\ miRNAs/Alderette-Goran\ EV\ Pilot/Datasets\ for\ Analysis/MOM155\ Removed/All\ RNAs\ Detected/Raw\ Counts/210914_Alderette_Goran_exceRpt_miRNAs_TMMnormalized_RawCounts_MOM155removed.csv")

# the first row and last 2 rows are not miR, the rest are
# 1219 - 3 = 1216

#read in the meta data
meta <- read.csv("/Volumes/IPHY/ADORLab/HEI Study/Master Datasets/Metadata/long/mothersMilk_metadata_timepointsAsRows_updated101921_Temporary24mDiet.csv")

#updated HMO code here: Z:
newHMO <- read.csv("/Volumes/IPHY/ADORLab/Lab\ Projects/Mothers\ Milk/HMO\ Master\ Data/Final_HMO.csv") 
#you'll need to merge that in after removing the old HMO variables from the meta data.

# Removing old and adding new HMOs
oldHMOs <- c("X2.FL..nmol.mL.", "X3FL..nmol.mL.", "LNnT..nmol.mL.",
"X3.SL..nmol.mL.", "DFLac..nmol.mL.", "X6.SL..nmol.mL.", "LNT..nmol.mL.",
"LNFP.I..nmol.mL.", "LNFP.II..nmol.mL.", "LNFP.III..nmol.mL.",
"LSTb..nmol.mL.", "LSTc..nmol.mL.", "DFLNT..nmol.mL.",
"LNH..nmol.mL.", "DSLNT..nmol.mL.", "FLNH..nmol.mL.",
"DFLNH..nmol.mL.", "FDSLNH..nmol.mL.", "DSLNH..nmol.mL.",
"SUM..nmol.mL.")
newHOMs <- c("x2FL_nmol_ml","x3FL_nmol_ml","LNnT_nmol_ml","x3SL_nmol_ml","DFLac_nmol_ml",
             "x6SL_nmol_ml","LNT_nmol_ml","LNFP_I_nmol_ml","LNFP_II_nmol_ml","LNFP_III_nmol_ml",
             "LSTb_nmol_ml","LSTc_nmol_ml","DFLNT_nmol_ml","LNH_nmol_ml","DSLNT_nmol_ml",
             "FLNH_nmol_ml","DFLNH_nmol_ml","FDSLNH_nmol_ml","DSLNH_nmol_ml","SUM_nmol_ml",
             "Sia_nmol_ml","Fuc_nmol_ml","x2FL_ug_ml","x3FL_ug_ml","LNnT_ug_ml",
             "x3SL_ug_ml", "DFLac_ug_ml", "x6SL_ug_ml", "LNT_ug_ml", "LNFP_I_ug_ml",
             "LNFP_II_percent","LNFP_III_percent","LSTb_percent","LSTc_percent","DFLNT_percent",
             "LNH_percent", "DSLNT_percent", "FLNH_percent", "DFLNH_percent", "FDSLNH_percent", 
             "DSLNH_percent")

meta_noold <- meta[,-which(names(meta) %in% oldHMOs)]
meta_new <- merge(meta, newHMO[,c("dyad_id", newHMOs)], by = "dyad_id", all.x = TRUE)

# only keep baseline variables for this analysis
table(meta$timepoint)
meta <- meta[meta$timepoint == 1,]

# CALCULATE DERIVED VARIABLES --------------------------------------------------

#breastfeedingcat####
#low
meta$breastfeedingcat[meta$breastmilk_per_day == 0] <- "None"
meta$breastfeedingcat[meta$breastmilk_per_day > 0 &
                        meta$breastmilk_per_day <= 3] <- "Low"
meta$breastfeedingcat[meta$breastmilk_per_day > 3 & 
                        meta$breastmilk_per_day < 8] <- "Med"
meta$breastfeedingcat[meta$breastmilk_per_day >= 8] <- "High"


#factor
meta$breastfeedingcat <- factor(meta$breastfeedingcat, levels = c("High", "Med",
                                                                  "Low", "None"))

summary(meta$breastfeedingcat)

#breastfeedings_continuous####
meta$breastfeedings_continuous <- 
  ifelse(meta$breastmilk_per_day == 0,0,
         ifelse(meta$breastmilk_per_day == 1,1,
                ifelse(meta$breastmilk_per_day > 1 &
                         meta$breastmilk_per_day < 9, 
                       meta$breastmilk_per_day-1, 8)))

hist(meta$breastfeedings_continuous)

#Diversity####
summary(meta$Diversity)
hist(meta$Diversity)

#Evenness####
summary(meta$Evenness)
hist(meta$Evenness)

#Sia####
summary(meta$Sia..nmol.mL.)
hist(meta$Sia..nmol.mL.)
#rename to Sia for simplicity
meta$Sia <- meta$Sia..nmol.mL.

#Fuc####
summary(meta$Fuc..nmol.mL.)
hist(meta$Fuc..nmol.mL.)

#rename to Fuc for simplicity
meta$Fuc <- meta$Fuc..nmol.mL.

#Secretor####
table(meta$Secretor)

# try to figure out what the ? status is
print(meta$X2.FL..nmol.mL.[meta$Secretor == "?"])

# see what the mean of 2'FL is for secretors vs. non
summary(meta$X2.FL..nmol.mL.[meta$Secretor == 1])
mean(meta$X2.FL..nmol.mL.[meta$Secretor == 0])

#clean up Secretor
meta$Secretor <-
  ifelse(meta$Secretor == "1", "Yes",
         ifelse(meta$Secretor == "0", "No", NA))

#factor
meta$Secretor <- factor(meta$Secretor)
summary(meta$Secretor)

#breast_milk_time_hrs####
summary(meta$breast_milk_time)
#time of collection is in military time
#convert to hours after midnight to make into continuous variable
meta$hours <- as.numeric(as.character(sub(":.*", "", 
                                          meta$breast_milk_time)))
meta$min <- as.numeric(as.character(sub(".*:", "", 
                                        meta$breast_milk_time)))

meta$breast_milk_time_hrs <- meta$hours + (meta$min/60)
head(meta$breast_milk_time_hrs, n = 5)
head(as.character(meta$breast_milk_time))

#mom_age_at_birth####
summary(meta$mom_age_at_birth)
hist(meta$mom_age_at_birth)

#SES####
sum(is.na(meta$SES_index_final))

#set missing values to the median SES
meta$SES_index_final[is.na(meta$SES_index_final)] <- 
  median(meta$SES_index_final, na.rm = T)
meta$SES <- meta$SES_index_final

summary(meta$SES)
hist(meta$SES)

#prepreg_bmi_kgm2####
summary(meta$prepreg_bmi_kgm2)
hist(meta$prepreg_bmi_kgm2)

#mom_BMI####
summary(meta$mom_BMI)
hist(meta$mom_BMI)

#M_fib_Mom####
summary(meta$m_fib_Mom)
hist(meta$m_fib_Mom)

#m_fat_Mom####
summary(meta$m_fat_Mom)
hist(meta$m_fat_Mom)

#m_pro_Mom####
summary(meta$m_pro_Mom)
hist(meta$m_pro_Mom)

#m_cho_Mom####
summary(meta$m_cho_Mom)
hist(meta$m_cho_Mom)

#m_fruc_Mom####
summary(meta$m_fruc_Mom)
hist(meta$m_fruc_Mom)

#m_asug_Mom####
summary(meta$m_asug_Mom)
hist(meta$m_asug_Mom)

#healthy_eating_ind####
hist(meta$HEI2015_TOTAL_SCORE_Mom)
meta$healthy_eating_ind <- meta$HEI2015_TOTAL_SCORE_Mom

#diet_inflam_ind####
hist(meta$DII_Mom)
meta$diet_inflam_ind <- meta$DII_Mom

#med_diet_score####
hist(meta$MDS_Mom)
meta$med_diet_score <- meta$MDS_Mom

# HMOs ####
hist(meta$X2.FL..nmol.mL.)
hist(meta$X3FL..nmol.mL.)
hist(meta$LNnT..nmol.mL.)
hist(meta$X3.SL..nmol.mL.)
hist(meta$DFLac..nmol.mL.)
hist(meta$X6.SL..nmol.mL.)
hist(meta$LNT..nmol.mL.)
hist(meta$LNFP.I..nmol.mL.)
hist(meta$LNFP.II..nmol.mL.)
hist(meta$LNFP.III..nmol.mL.)
hist(meta$LSTb..nmol.mL.)
hist(meta$LSTc..nmol.mL.)
hist(meta$DFLNT..nmol.mL.)
hist(meta$LNH..nmol.mL.)
hist(meta$DSLNT..nmol.mL.)
hist(meta$FLNH..nmol.mL.)
hist(meta$DFLNH..nmol.mL.)
hist(meta$FDSLNH..nmol.mL.)
hist(meta$DSLNH..nmol.mL.)
hist(meta$SUM..nmol.mL.)

#baby_gender_cat ####
meta$baby_gender_cat <- factor(ifelse(meta$baby_gender %in% 1, 
                                      "Female","Male"))
summary(meta$baby_gender_cat)

#gestational_age_cat####
#(on-time, early [>2 weeks before due date], late [>2 weeks after due date]) 
meta$gestational_age_cat <- 
  factor(ifelse(meta$gestational_age_category %in% c("<38","38-40"),
                "Early", 
                ifelse(meta$gestational_age_category %in% c("40-42",">42"),
                       "Late","On Time")),
         levels = c("On Time","Late","Early"))

summary(meta$gestational_age_cat)

#rapidGrowth####
meta$rapidGrowth <- factor(meta$rapidGrowth)
summary(meta$rapidGrowth)

# predicted fat mass####
#predicted fat mass (have to do on complete cases so we don't get negative vals)
meta$pred_fat_mass <- -0.61556255 + (4.09644221*I((meta$inf_weight_kg/10)^3)) - 
  (6.87699828*I((meta$inf_weight_kg/10)^3*log((meta$inf_weight_kg/10)))) - 
  (3.58838763*I((meta$inf_length_cm/100)^3)) - 
  (14.70363829*I((meta$inf_length_cm/100)^3*log((meta$inf_length_cm/100)))) - 
  (0.03962223*I((meta$age_in_days/100)^3)) + 
  (0.02136099*I((meta$age_in_days/100)^3 * log((meta$age_in_days/100)))) - 
  (0.18975680*meta$baby_gender)

hist(meta$pred_fat_mass)
print(meta$merge_id_dyad[which(meta$pred_fat_mass < 0)])
print(meta$merge_id_dyad[which(meta$pred_fat_mass > 50)])

# don't worry about the crazy measures - they are at 24 and 36 months and we
# won't use those

# TSF ####
#convert from string to numeric
meta$skinf_midthigh_mm <- as.numeric(as.character(meta$skinf_midthigh_mm))
sum(is.na(meta$skinf_midthigh_mm))
meta$skinf_tricep_mm <- as.numeric(as.character(meta$skinf_tricep_mm))
meta$skinf_subscap_mm <- as.numeric(as.character(meta$skinf_subscap_mm))
meta$skinf_supra_mm <- as.numeric(as.character(meta$skinf_supra_mm))

meta$TSF <- meta$skinf_midthigh_mm + meta$skinf_tricep_mm +
  meta$skinf_subscap_mm + meta$skinf_supra_mm

hist(meta$TSF)
sum(is.na(meta$TSF))

#Calculate CTSF

meta$CTSF <- (meta$skinf_subscap_mm + meta$skinf_supra_mm)/meta$TSF

hist(meta$CTSF)

#mom_abx####
table(meta$mom_antibiotics)

#baby_abx####
table(meta$baby_antibiotics)

meta$baby_abx <- meta$baby_antibiotics

meta$baby_abx <- na_if(meta$baby_abx, 4)
sum(is.na(meta$baby_abx))

meta$baby_abx <- if_else(meta$baby_abx > 1, 1, 0)

meta$baby_abx <- factor(meta$baby_abx, labels = c("No", "Yes"))

summary(meta$baby_abx)

#mode_of_delivery_cat####
meta$mode_of_delivery_cat <- 
  factor(ifelse(meta$mode_of_delivery %in% 1, 
                "Vaginal","C-Section"),
         levels = c("Vaginal","C-Section"))
summary(meta$mode_of_delivery_cat)

# season ####
head(meta$date_of_visit)

# convert to a date 
meta$visit_date <- mdy(meta$date_of_visit)
summary(meta$visit_date)

meta$season <- ifelse(month(meta$visit_date) %in% c(10, 11, 12, 1, 2, 3), 
                      "Cold", "Warm")

meta$season <- as.factor(meta$season)

# formula_reg####

summary(meta$formula_reg)

meta$formula_reg <- factor(meta$formula_reg, labels = c("Yes", "No"))
summary(meta$formula_reg)

# form_vs_breast ####
summary(meta$form_vs_breast)
hist(meta$form_vs_breast)

# formula_continuous ####
meta$formula_continuous <- 
  ifelse(meta$formula_per_day == 0,0,
         ifelse(meta$formula_per_day == 1,1,
                ifelse(meta$formula_per_day > 1 &
                         meta$formula_per_day < 9, 
                       meta$formula_per_day-1, 8)))
summary(meta$formula_continuous)

# pred_bf ####
meta$pred_bf <- factor(if_else(meta$formula_reg == "Yes", "No", "Yes"))
summary(meta$pred_bf)

# formulacat ####
meta$formulacat[meta$formula_continuous == 0] <- "None"
meta$formulacat[meta$formula_continuous > 0 & 
                  meta$formula_continuous <= 3] <- "Low"
meta$formulacat[meta$formula_continuous > 3 & meta$formula_continuous < 8] <- "Med"
meta$formulacat[meta$formula_continuous >= 8] <- "High"

meta$formulacat <- factor(meta$formulacat, levels = c("High", "Med", "Low",
                                                      "None"))

summary(meta$formulacat)

#subset on only the variables we need
meta_trim <- dplyr::select(meta, c(dyad_id, merge_id_dyad, breastfeedingcat,
                                   breastfeedings_continuous,
                                   Diversity, Evenness, Sia, Fuc,
                                   Secretor, breast_milk_time_hrs,
                                   mom_age_at_birth, SES,
                                   prepreg_bmi_kgm2, mom_BMI,
                                   healthy_eating_ind, 
                                   diet_inflam_ind, m_ener_Mom,
                                   med_diet_score, baby_gender_cat,
                                   gestational_age_cat, rapidGrowth, 
                                   baby_abx, mode_of_delivery_cat,
                                   age_in_days, 
                                   m_fat_Mom, m_pro_Mom, m_cho_Mom,
                                   m_fruc_Mom, m_asug_Mom, inf_weight_kg,
                                   inf_length_cm, zwfl, zwei, circ_umb_cm,
                                   TSF, CTSF, m_fib_Mom, X2.FL..nmol.mL.,
                                   X3FL..nmol.mL., LNnT..nmol.mL.,
                                   X3.SL..nmol.mL., DFLac..nmol.mL.,
                                   X6.SL..nmol.mL., LNT..nmol.mL.,
                                   LNT..nmol.mL., LNFP.I..nmol.mL.,
                                   LNFP.II..nmol.mL., LNFP.III..nmol.mL.,
                                   LSTb..nmol.mL., LSTc..nmol.mL., 
                                   DFLNT..nmol.mL., LNH..nmol.mL.,
                                   DSLNT..nmol.mL., FLNH..nmol.mL.,
                                   DFLNH..nmol.mL., FDSLNH..nmol.mL.,
                                   DSLNH..nmol.mL., SUM..nmol.mL.,
                                   mom_antibiotics, season, pred_fat_mass,
                                   formula_reg, form_vs_breast, 
                                   formula_continuous, formulacat, pred_bf))

# save the meta data
write.csv(meta_trim, file = paste0(rda_out, "meta_clean_bl.csv"))
#meta_clean_bl.csv####

# CLEAN MIRNA DATA -------------------------------------------------------------
# Allison suggests looking at the percent RNA, supernatant 
# volume, and EV isolation date
hist(miRNA_QC$prop_rRNA)
hist(miRNA_QC$Vol_Supernatant)
summary(miRNA_QC$Date_Evs)
summary(miRNA_QC$Low_TranscriptomeGenomeRatio)

# convert Date_Evs to factor
miRNA_QC$Date_Evs <- as.factor(miRNA_QC$Date_Evs)

#subset on the variables that we will need to adjust for
miRNA_QC_sub <- dplyr::select(miRNA_QC, c("Sample_ID", "prop_rRNA",
                                          "Vol_Supernatant",
                                          "Date_Evs", "prop_unmapped",
                                          "Low_TranscriptomeGenomeRatio"))

#update the miRNA_QC_sub row names that don't have MOM in front of them
miRNA_QC_sub$Sample_ID <- as.character(miRNA_QC_sub$Sample_ID)

miRNA_QC_sub$Sample_ID[miRNA_QC_sub$Sample_ID == "166"] <- "MOM.166"
miRNA_QC_sub$Sample_ID[miRNA_QC_sub$Sample_ID == "168"] <- "MOM.168"

# updates to miRNA_cpm ####
#transpose
miRNA_t <- data.frame(t(miRNA_cpm))

#update the column names
miRNA_names <- t(miRNA_cpm$X)

colnames(miRNA_t) <- miRNA_names

miRNA_t <- miRNA_t[-1,]

#update X166 and X168 to be the same format as the other variables
row.names(miRNA_t)[row.names(miRNA_t) == "X166"] <- "MOM.166"
row.names(miRNA_t)[row.names(miRNA_t) == "X168"] <- "MOM.168"

miRNA_t$Sample_ID <- row.names(miRNA_t)

#join in the QC variables to miRNA_t
miRNA_t <- left_join(miRNA_t, miRNA_QC_sub, by = "Sample_ID")

#make a dyad_id variable
miRNA_t$dyad_id <- substr((miRNA_t$Sample_ID), 5, 7)

#drop the leading zeros from dyad id
miRNA_t$dyad_id <- as.numeric(str_remove(miRNA_t$dyad_id, "^0+"))

#make the dyad_id the row name
row.names(miRNA_t) <- miRNA_t$dyad_id

#make all variables in miRNA_t numeric
for(i in 1:210){ #the first 210 observations are the miRNA - the rest are meta
  miRNA_t[,i] <- as.numeric(as.character(miRNA_t[,i]))
}

#save the miRNA_t file
write.csv(miRNA_t, file = paste0(rda_out, "miRNA_cpm.csv"))
#miRNA_cpm.csv####

# ------------------------------------
# updates to miRNA_counts####

#update the first column name
colnames(miRNA_counts)[1] <- "Sample_ID"

#update X166 and X168 to be the same format as the other variables
miRNA_counts$Sample_ID <- as.character(miRNA_counts$Sample_ID)
miRNA_counts$Sample_ID[miRNA_counts$Sample_ID == "166"] <- "MOM.166"
miRNA_counts$Sample_ID[miRNA_counts$Sample_ID == "168"] <- "MOM.168"

#join in the QC variables to miRNA_t
miRNA_counts <- left_join(miRNA_counts, miRNA_QC, by = "Sample_ID")

#make a dyad_id variable
miRNA_counts$dyad_id <- substr((miRNA_counts$Sample_ID), 5, 7)

#drop the leading zeros from dyad id
miRNA_counts$dyad_id <- as.numeric(str_remove(miRNA_counts$dyad_id, 
                                              "^0+"))

#make the dyad_id the row name
row.names(miRNA_counts) <- miRNA_counts$dyad_id

#make all variables in miRNA_t numeric
for(i in 2:211){ #the first 210 observations are the miRNA - the rest are meta
  miRNA_counts[,i] <- as.numeric(as.character(miRNA_counts[,i]))
}

#save the clean miRNA counts file
write.csv(miRNA_counts, file = paste0(rda_out, "miRNA_counts.csv"))
#miRNA_counts.csv####

# miRNA SUMMARY STATS ----------------------------------------------------------

# get the number of miRNA present in all samples
perc.zero.raw <- data.frame(matrix(ncol = 210, nrow = 1))
colnames(perc.zero.raw) <- colnames(miRNA_counts[,2:211])
perc.zero.raw <- colSums(miRNA_counts[,2:211] == 0)/110
sum(perc.zero.raw == 0) #80
hist(perc.zero.raw, main = "Percent 0's", xlab = "")

# get the total number of reads
sum(miRNA_counts$LibSize) #171991022
# get average reads per sample #1563555
sum(miRNA_counts$LibSize)/nrow(miRNA_counts)

# sort by colSums to get the top 5 most abundant
miRNA_counts_abund <- miRNA_counts[,2:211]
miRNA_counts_abund <- miRNA_counts_abund[,order(colSums(miRNA_counts_abund))]

# TABLE 1 ----------------------------------------------------------------------
# make a list of vars to summarize
summary(meta_trim)
t1_vars <- c("dyad_id", "mom_age_at_birth", "SES", "mom_BMI",  "baby_gender_cat",  
             "mode_of_delivery_cat", "age_in_days", "breast_milk_time_hrs", 
             "pred_bf", "formula_continuous" ,"breastfeedings_continuous", 
             "Diversity", "Sia",
             "Fuc", "Secretor", "X2.FL..nmol.mL.", "X3FL..nmol.mL.", "LNnT..nmol.mL.",
             "X3.SL..nmol.mL.", "DFLac..nmol.mL.", "X6.SL..nmol.mL.", "LNT..nmol.mL.",
             "LNFP.I..nmol.mL.", "LNFP.II..nmol.mL.", "LNFP.III..nmol.mL.",
             "LSTb..nmol.mL.", "LSTc..nmol.mL.", "DFLNT..nmol.mL.",
             "LNH..nmol.mL.", "DSLNT..nmol.mL.", "FLNH..nmol.mL.",
             "DFLNH..nmol.mL.", "FDSLNH..nmol.mL.", "DSLNH..nmol.mL.",
             "SUM..nmol.mL.")

t1_vars <- meta_trim[,colnames(meta_trim) %in% t1_vars]

#subset on those with miRNA data
#t1_vars <- t1_vars[t1_vars$dyad_id %in% row.names(miRNA_t),]
t1_vars <- t1_vars[t1_vars$dyad_id %in% row.names(miRNA_counts),]

# summarize the missingness
missing_plot(t1_vars)

# missing secretor status for 2 individuals

# output continuous summary stats
stargazer(t1_vars, type = "text")

# now summarize the categorical variables
summary(t1_vars$baby_gender_cat)
summary(t1_vars$mode_of_delivery_cat)
summary(t1_vars$pred_bf)
summary(t1_vars$Secretor)

hist(t1_vars$breastfeedings_continuous)

# get the number of feedings per day among those who are predominantly breast fed
mean(t1_vars$breastfeedings_continuous[which(t1_vars$pred_bf == "Yes")])
sd(t1_vars$breastfeedings_continuous[which(t1_vars$pred_bf == "Yes")])

# CORRELATION ------------------------------------------------------------------
cor.test(t1_vars$formula_continuous, t1_vars$breastfeedings_continuous)
t1_vars$formula_reg <- ifelse(t1_vars$formula_reg == "Yes", 1, 0)
cor.test(t1_vars$formula_reg, t1_vars$breastfeedings_continuous)

# make a correlation plot for the HMOs 
corr_data <- dplyr::select(t1_vars, c("Sia", "Fuc", "SUM..nmol.mL.", 
"X2.FL..nmol.mL.", "X3FL..nmol.mL.", "X3.SL..nmol.mL.", 
"X6.SL..nmol.mL.", "DFLac..nmol.mL.", "DFLNH..nmol.mL.", 
"DFLNT..nmol.mL.", "DSLNH..nmol.mL.", "DSLNT..nmol.mL.",
"FDSLNH..nmol.mL.", "FLNH..nmol.mL.", "LNnT..nmol.mL.", 
"LNT..nmol.mL.", "LNFP.I..nmol.mL.", "LNFP.II..nmol.mL.", 
"LNFP.III..nmol.mL.", "LNH..nmol.mL.", "LSTb..nmol.mL.", 
"LSTc..nmol.mL."))

# make labels for the correlation plot
colnames(corr_data) <- c("Sialyated HMOs", "Fucosylated HMOs", "Sum of HMOs",
          "2'-fucosyllactose", "3'-fucosyllactose", "3'-sialyllactose",
          "6'-sialyllactose", "DFLac", "DFLNH", "DFLNT", "DSLNH", "DSLNT",
          "FDSLNH", "FLNH", "Lacto-N-neotetraose", "Lacto-N-tetraose",
          "LNFP I", "LNFP II", "LNFP III", "LNH", "LSTb", "LSTc")

M = cor(corr_data)
testRes = cor.mtest(corr_data, conf.level = 0.95)

png(filename = paste0(figs_out, "Figure_1.png"), 
    units = "in",
    width = 10,
    height = 8,
    pointsize = 12,
    res = 1200)
corrplot(M, method = "color", diag = T, type = "upper", p.mat = testRes$p,
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, insig = 'label_sig', 
         pch.col = "black", tl.col='black')
dev.off()

# Figure_1.png ####
