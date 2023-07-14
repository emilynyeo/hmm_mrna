# HEADER -----------------------------------------------------------------------
#
# TITLE:   3-cluster-analysis.R
#
# PURPOSE: make clusters based on miRNA expression, analyze whether there are
#          differences in HMO concentration between clusters
#
# DATE:    March 23, 2022
#
# SET UP -----------------------------------------------------------------------

#clear the workspace
rm(list = ls())

#turn off scientific notation
options(scipen = 100)

#load libraries
library(plyr); library(tidyverse); library(factoextra); library(xlsx)
library(purrr); library(stringr); library(lme4); library(dendextend)
library(lmerTest); library(corrplot); library(lubridate)
library(MASS); library(pscl)

#set the input folder
data_in <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/rda"

#set the output folder
figs_out <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/figs"

#read in the clean miRNA_cpm data
miRNA_cpm <- read.csv(file = paste0(data_in, "miRNA_counts.csv"))

#read in the clean meta data
meta <- read.csv(file = paste0(data_in, "meta_clean_bl.csv"))

# PREPARE THE DATA SET ---------------------------------------------------------
# drop unnecessary columns
rownames(miRNA_cpm) <- paste0("X", miRNA_cpm$dyad_id)

#drop meta data
miRNA_cpm_trim <- dplyr::select(miRNA_cpm, -c("X", "Sample_ID", "prop_rRNA",
                                              "Vol_Supernatant", "Date_Evs",
                                              "prop_unmapped", 
                                              "Low_TranscriptomeGenomeRatio",
                                              "dyad_id"))

# make a second miRNA cpm object with only secretors
table(meta$Secretor)
secretors <- meta$dyad_id[which(meta$Secretor == "Yes")]

# only get miRNA of secretors to test whether clustering is better
miRNA_sec <- miRNA_cpm[which(miRNA_cpm$X %in% secretors),]
#96 of 110 were secretors

# drop the extra variables
miRNA_sec_trim <- miRNA_sec[,2:211]

# scale 
miRNA_cpm_scaled <- data.frame(scale(miRNA_cpm_trim))
miRNA_sec_trim <- data.frame(scale(miRNA_sec_trim))

# now calculate the distance matrix
dist_mat <- dist(miRNA_cpm_scaled, method = 'euclidean')
dist_mat_sec <- dist(miRNA_sec_trim, method = "euclidean")

# cluster
hclust_avg <- hclust(dist_mat, method = "complete")
hclust_sec <- hclust(dist_mat_sec, method = "complete")

# make a plot
dend <- as.dendrogram(hclust_avg)
dend_sec <- as.dendrogram(hclust_sec)

dend %>% plot

dend_sec %>% plot 

#change the color of the branches to correspond with the clusters
dend %>% set("branches_k_color", k = 9) %>%
  set("labels_cex", 0) %>%
  plot()

plot(hclust_avg, cex = 0.4)

plot(hclust_sec, cex = 0.4)

# use the "Elbow" method to estimate the appropriate number of clusters
set.seed(78)
fviz_nbclust(miRNA_cpm_scaled, kmeans, method = "wss", k.max = 10)

fviz_nbclust(miRNA_sec_trim, kmeans, method = "wss", k.max = 10)
# 3 or 4 clusters,  based on this

# we should use 9 clusters based on the figure, but that results in several
# groups having only 1 person in them - try 3 and see how that looks
# with three, there is still only one individual in the third cluster (X166)
# try 2

hclust <- cutree(hclust_avg, 6)

hclust_sec2 <- cutree(hclust_sec, 4)

# get some info about the clustering
table(hclust)
table(hclust_sec2)

# the clustering is better if we only consider secretors - do it both ways

# use the original hclust for now, comparing clusters 1 and 3

#make a list of individuals in cluster 1 and cluster 3
cluster_1 <- hclust[hclust == 1]
cluster_3 <- hclust[hclust == 3]

# add cluster as a variable in meta
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_1)] <- 1
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_3)] <- 3

summary(meta$cluster)

hc <- meta[which(!is.na(meta$cluster)),]

# COMPARE CLUSTERS 1 AND 3 -----------------------------------------------------

# make a list of vars to include
vars <- c("mom_age_at_birth", "SES", "baby_gender_cat", "mode_of_delivery_cat", 
          "age_in_days", "pred_bf", "breastfeedings_continuous", "mom_BMI",
          "Secretor", "Diversity", "Sia", "Fuc", "SUM..nmol.mL.", 
          "X2.FL..nmol.mL.", "X3FL..nmol.mL.", "X3.SL..nmol.mL.", 
          "X6.SL..nmol.mL.", "DFLac..nmol.mL.", "DFLNH..nmol.mL.", 
          "DFLNT..nmol.mL.", "DSLNH..nmol.mL.", "DSLNT..nmol.mL.",
          "FDSLNH..nmol.mL.", "FLNH..nmol.mL.", "LNnT..nmol.mL.", 
          "LNT..nmol.mL.", "LNFP.I..nmol.mL.", "LNFP.II..nmol.mL.", 
          "LNFP.III..nmol.mL.", "LNH..nmol.mL.", "LSTb..nmol.mL.", 
          "LSTc..nmol.mL.", "breast_milk_time_hrs")
vars <- unlist(vars)

# make separate data sets for clusters 1 and 3
cluster_1 <- hc[which(hc$cluster == 1),]
cluster_3 <- hc[which(hc$cluster == 3),]

# make an output table
cluster_out <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(cluster_out) <- c("Var", "Cluster_1", "Cluster_3", "p")

# loop over vars, run t-tests or chi-squared tests for each
for(thisVar in vars){
  # set the value for the next row
  i <- nrow(cluster_out) + 1
  # categorical variables
  if(thisVar %in% c("baby_gender_cat", "mode_of_delivery_cat", "pred_bf",
                     "Secretor")) { 
    temp <- chisq.test(hc[,thisVar], hc$cluster)
    cluster_out[i, "Var"] <- thisVar
    cluster_out[i, "p"] <- temp$p.value
  }
  # continuous variables
  else{ 
    temp <- t.test(cluster_1[,thisVar], cluster_3[,thisVar])
    # add the results to the output table
    cluster_out[i, "Var"] <- thisVar
    cluster_out[i, "Cluster_1"] <- paste0(round(temp$estimate[1], digits = 2), 
                                          " ± ", 
                                          round(sd(cluster_1[,thisVar]), 
                                                   digits = 2))
    cluster_out[i, "Cluster_3"] <- paste0(round(temp$estimate[2], digits = 2), 
                                          " ± ",
                                          round(sd(cluster_3[,thisVar]), digits = 2))
    cluster_out[i, "p"] <- temp$p.value
  }
}

# correct for multiple testing
cluster_out$pfdr <- p.adjust(cluster_out$p)

# get the n for the categorical variables in each cluster
table(cluster_1$baby_gender_cat)
table(cluster_3$baby_gender_cat)

table(cluster_1$mode_of_delivery_cat)
table(cluster_3$mode_of_delivery_cat)

table(cluster_1$pred_bf)
table(cluster_3$pred_bf)

table(cluster_1$Secretor)
table(cluster_3$Secretor)


# SENSITIVITY ANALYSIS - SECRETORS ---------------------------------------------

#make a list of individuals in cluster 1 and cluster 3
cluster_1 <- hclust_sec2[hclust_sec2 == 1]
cluster_2 <- hclust_sec2[hclust_sec2 == 2]
cluster_3 <- hclust_sec2[hclust_sec2 == 3]

# add cluster as a variable in meta
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_1)] <- 1
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_2)] <- 2
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_3)] <- 3

summary(meta$cluster)

hc <- meta[which(!is.na(meta$cluster)),]

# remove secretor from vars
vars <- c("mom_age_at_birth", "SES", "baby_gender_cat", "mode_of_delivery_cat", 
          "age_in_days", "pred_bf", "breastfeedings_continuous", "mom_BMI",
          "Diversity", "Sia", "Fuc", "SUM..nmol.mL.", 
          "X2.FL..nmol.mL.", "X3FL..nmol.mL.", "X3.SL..nmol.mL.", 
          "X6.SL..nmol.mL.", "DFLac..nmol.mL.", "DFLNH..nmol.mL.", 
          "DFLNT..nmol.mL.", "DSLNH..nmol.mL.", "DSLNT..nmol.mL.",
          "FDSLNH..nmol.mL.", "FLNH..nmol.mL.", "LNnT..nmol.mL.", 
          "LNT..nmol.mL.", "LNFP.I..nmol.mL.", "LNFP.II..nmol.mL.", 
          "LNFP.III..nmol.mL.", "LNH..nmol.mL.", "LSTb..nmol.mL.", 
          "LSTc..nmol.mL.", "breast_milk_time_hrs")
vars <- unlist(vars)

# make separate data sets for clusters 1 and 3
cluster_1 <- hc[which(hc$cluster == 1),]
cluster_2 <- hc[which(hc$cluster == 2),]
cluster_3 <- hc[which(hc$cluster == 3),]

hc$cluster <- factor(hc$cluster)

# make an output table
cluster_sec <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(cluster_sec) <- c("Var", "Mean_SD_C1", "Mean_SD_C2", "Mean_SD_C3",
                           "p2", "p3")

# loop over vars, run t-tests or chi-squared tests for each
for(thisVar in vars){
  # set the value for the next row
  i <- nrow(cluster_sec) + 1
  # categorical variables
  if(thisVar %in% c("baby_gender_cat", "mode_of_delivery_cat", "pred_bf")) { 
    temp <- chisq.test(hc[,thisVar], hc$cluster)
    cluster_sec[i, "Var"] <- thisVar
    cluster_sec[i, "p"] <- temp$p.value
  }
  # continuous variables
  else{ 
    # add the results to the output table
    cluster_sec[i, "Var"] <- thisVar
    cluster_sec[i, "Mean_SD_C1"] <- paste(round(mean(cluster_1[,thisVar]), digits = 4), 
                                          "±",
                                          round(sd(cluster_1[,thisVar]), digits = 4), sep = " ")
    cluster_sec[i, "Mean_SD_C2"] <- paste(round(mean(cluster_2[,thisVar]), digits = 4), 
                                          "±",
                                          round(sd(cluster_2[,thisVar]), digits = 4), sep = " ")
    cluster_sec[i, "Mean_SD_C3"] <- paste(round(mean(cluster_3[,thisVar]), digits = 4), 
                                          "±",
                                          round(sd(cluster_3[,thisVar]), digits = 4), sep = " ")
    cluster_sec[i, "p2"] <- round(t.test(cluster_1[,thisVar], 
                                         cluster_2[,thisVar])$p.value, digits = 4)
    cluster_sec[i, "p3"] <- round(t.test(cluster_1[,thisVar], 
                                         cluster_3[,thisVar])$p.value, digits = 4)
  }
}

# correct for multiple testing
cluster_sec$pfdr <- p.adjust(cluster_sec$p)
cluster_sec$pfdr2 <- p.adjust(cluster_sec$p2)
cluster_sec$pfdr3 <- p.adjust(cluster_sec$p3)


