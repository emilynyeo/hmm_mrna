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
pacman::p_load(knitr, tidyverse, magrittr, lme4, lmerTest, GGally, corrplot, 
               Hmisc, kableExtra, dplyr, plyr, janitor, lubridate, survminer, 
               ggplot2, here, readr, tableone, officer, flextable,finalfit,
               purrr, stringr, lme4, corrplot, pscl, stargazer, MASS, lmerTest,
               readxl, factoextra, ggbiplot, dendextend)

#set the input folder
#data_in <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/rda"

#set the output folder
#figs_out <- "/Volumes/IPHY/ADORLab/__Users/emye7956/MM/HMO-miRNA/1-data-cleaning/figs"

#read in the clean miRNA_cpm data
#miRNA_cpm <- read.csv(file = paste0(data_in, "miRNA_counts.csv"))
miRNA_cpm <- read.csv("input/miRNA_counts_EY.csv")

#read in the clean meta data
#meta <- read.csv(file = paste0(data_in, "meta_clean_bl.csv"))
meta <- read.csv("input/meta_clean_EY.csv")

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
#No Yes 
#24 183 
secretors <- meta$dyad_id[which(meta$Secretor == "Yes")]

# only get miRNA of secretors to test whether clustering is better
miRNA_sec <- miRNA_cpm[which(miRNA_cpm$X %in% secretors),]
#96 of 110 were secretors

# drop the extra variables
miRNA_sec_trim_ID <- miRNA_sec[,2:211]
#remove sample names 
miRNA_sec_trim_ID <- miRNA_sec_trim_ID[, -1]
#make numeric
miRNA_cpm_trim <- sapply(miRNA_cpm_trim, as.numeric)
miRNA_sec_trim <- sapply(miRNA_sec_trim_ID, as.numeric)

# scale 
miRNA_cpm_scaled <- data.frame(scale(miRNA_cpm_trim))
miRNA_sec_trim <- data.frame(scale(miRNA_sec_trim))

# now calculate the distance matrix
dist_mat <- dist(miRNA_cpm_scaled, method = 'euclidean')
dist_mat_sec <- dist(miRNA_sec_trim, method = "euclidean")

# cluster
hclust_avg <- hclust(dist_mat, method = "complete")
hclust_sec <- hclust(dist_mat_sec, method = "complete")
#revise methods complete 

# make a plot
dend <- as.dendrogram(hclust_avg)
dend_sec <- as.dendrogram(hclust_sec)

dend %>% plot

dend_sec %>% plot 

#change the color of the branches to correspond with the clusters

#CHANGE TO SEE WHAT THIS LOOKS LIKE CLUSTERING BY HMO
dend %>% set("branches_k_color", k = 9) %>%
  set("labels_cex", 0) %>%
  plot()

plot(hclust_avg, cex = 0.4)

plot(hclust_sec, cex = 0.4)

# use the "Elbow" method to estimate the appropriate number of clusters
set.seed(78)
# Remove columns with NA values
miRNA_cpm_scaled <- miRNA_cpm_scaled[, colSums(is.na(miRNA_cpm_scaled)) == 0]
fviz_nbclust(miRNA_cpm_scaled, kmeans, method = "wss", k.max = 10)
fviz_nbclust(miRNA_sec_trim, kmeans, method = "wss", k.max = 10)
# 3 or 4 clusters,  based on this

# we should use 9 clusters based on the figure, but that results in several
# groups having only 1 person in them - try 3 and see how that looks
# with three, there is still only one individual in the third cluster (X166)
# try 2

hclust <- cutree(hclust_avg, 6)
names(hclust) <- miRNA_cpm$dyad_id
hclust_sec2 <- cutree(hclust_sec, 4)
#names(hclust_sec2) <- miRNA_cpm$dyad_id

# get some info about the clustering
table(hclust)
#1  2  3  4  5  6 
#28 69  2  9  1  1 
table(hclust_sec2)
#1  2  3  4 
#16 77  2  1 

# the clustering is better if we only consider secretors - do it both ways

# use the original hclust for now, comparing clusters 1 and 3

#make a list of individuals in cluster 1 and cluster 3
cluster_1 <- hclust[hclust == 1]
cluster_2 <- hclust[hclust == 2]
cluster_3 <- hclust[hclust == 3]

# add cluster as a variable in meta
meta$cluster <- NA
meta$cluster[row.names(meta) %in% names(cluster_1)] <- 1
meta$cluster[row.names(meta) %in% names(cluster_2)] <- 2
meta$cluster[row.names(meta) %in% names(cluster_3)] <- 3

summary(meta$cluster)

hc <- meta[which(!is.na(meta$cluster)),]

# COMPARE CLUSTERS 1 AND 2 -----------------------------------------------------

# make a list of vars to include
vars <- c("mom_age_at_birth", "SES", "baby_gender_cat", "mode_of_delivery_cat", 
          "age_in_days", "pred_bf", "breastfeedings_continuous", "mom_BMI",
          "Secretor", "Diversity", "Sia", "Fuc", #"Secretor", 
          "x2FL_nmol_ml","x3FL_nmol_ml",
          "LNnT_nmol_ml",
          "x3SL_nmol_ml","DFLac_nmol_ml",
          "x6SL_nmol_ml","LNT_nmol_ml","LNFP_I_nmol_ml",
          "LNFP_II_nmol_ml","LNFP_III_nmol_ml",
          "LSTb_nmol_ml","LSTc_nmol_ml","DFLNT_nmol_ml",
          "LNH_nmol_ml",
          "DSLNT_nmol_ml","FLNH_nmol_ml","DFLNH_nmol_ml",
          "FDSLNH_nmol_ml",
          "DSLNH_nmol_ml","SUM_nmol_ml",
          "Sia_nmol_ml","Fuc_nmol_ml","x2FL_ug_ml",
          "x3FL_ug_ml","LNnT_ug_ml",
          "x3SL_ug_ml", "DFLac_ug_ml", "x6SL_ug_ml", 
          "LNT_ug_ml", "LNFP_I_ug_ml", "breast_milk_time_hrs")
#CHANGE OLD NAMES

vars <- unlist(vars)

# make separate data sets for clusters 1 and 3
cluster_1 <- hc[which(hc$cluster == 1),]
cluster_2 <- hc[which(hc$cluster == 2),]
cluster_3 <- hc[which(hc$cluster == 3),]

# make an output table
cluster_out <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(cluster_out) <- c("Var", "Cluster_1", "Cluster_2", "Cluster_3", "p")

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
    temp <- t.test(cluster_1[,thisVar], cluster_2[,thisVar])
    # add the results to the output table
    cluster_out[i, "Var"] <- thisVar
    cluster_out[i, "Cluster_1"] <- paste0(round(temp$estimate[1], digits = 2), 
                                          " ± ", 
                                          round(sd(cluster_1[,thisVar]), 
                                                   digits = 2))
    cluster_out[i, "Cluster_2"] <- paste0(round(temp$estimate[2], digits = 2), 
                                          " ± ",
                                          round(sd(cluster_2[,thisVar]), digits = 2))
    cluster_out[i, "Cluster_3"] <- paste0(round(temp$estimate[2], digits = 2), 
                                          " ± ",
                                          round(sd(cluster_3[,thisVar]), digits = 2))
    cluster_out[i, "p"] <- temp$p.value
  }
}

# correct for multiple testing
cluster_out$pfdr <- p.adjust(cluster_out$p) # no sig diff

# get the n for the categorical variables in each cluster
table(cluster_1$baby_gender_cat)
table(cluster_2$baby_gender_cat)
table(cluster_3$baby_gender_cat)

table(cluster_1$mode_of_delivery_cat)
table(cluster_2$mode_of_delivery_cat)
table(cluster_2$mode_of_delivery_cat)

table(cluster_1$pred_bf)
table(cluster_2$pred_bf)
table(cluster_3$pred_bf)

table(cluster_1$Secretor)
table(cluster_2$Secretor)
table(cluster_3$Secretor)

# SENSITIVITY ANALYSIS - SECRETORS ---------------------------------------------

#make a list of individuals in cluster 1 and cluster 3
cluster_1 <- hclust_sec2[hclust_sec2 == 1]
cluster_2 <- hclust_sec2[hclust_sec2 == 2]
cluster_3 <- hclust_sec2[hclust_sec2 == 3]
cluster_4 <- hclust_sec2[hclust_sec2 == 4]
cluster_5 <- hclust_sec2[hclust_sec2 == 5]

# add cluster as a variable in meta
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_1)] <- 1
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_2)] <- 2
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_3)] <- 3
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_4)] <- 4
meta$cluster[paste0("X", meta$dyad_id) %in% names(cluster_5)] <- 5

summary(meta$cluster)

hc <- meta[which(!is.na(meta$cluster)),]

# remove secretor from vars
vars <- c("mom_age_at_birth", "SES", "baby_gender_cat", "mode_of_delivery_cat", 
          "age_in_days", "pred_bf", "breastfeedings_continuous", "mom_BMI",
          "Diversity", "Sia", "Fuc", #"Secretor", 
          "x2FL_nmol_ml","x3FL_nmol_ml",
          "LNnT_nmol_ml",
          "x3SL_nmol_ml","DFLac_nmol_ml",
          "x6SL_nmol_ml","LNT_nmol_ml","LNFP_I_nmol_ml",
          "LNFP_II_nmol_ml","LNFP_III_nmol_ml",
          "LSTb_nmol_ml","LSTc_nmol_ml","DFLNT_nmol_ml",
          "LNH_nmol_ml",
          "DSLNT_nmol_ml","FLNH_nmol_ml","DFLNH_nmol_ml",
          "FDSLNH_nmol_ml",
          "DSLNH_nmol_ml","SUM_nmol_ml",
          "Sia_nmol_ml","Fuc_nmol_ml","x2FL_ug_ml",
          "x3FL_ug_ml","LNnT_ug_ml",
          "x3SL_ug_ml", "DFLac_ug_ml", "x6SL_ug_ml", 
          "LNT_ug_ml", "LNFP_I_ug_ml", "breast_milk_time_hrs")
vars <- unlist(vars)

# make separate data sets for clusters 1 and 3
cluster_1 <- hc[which(hc$cluster == 1),]
cluster_2 <- hc[which(hc$cluster == 2),]
cluster_3 <- hc[which(hc$cluster == 3),]
cluster_4 <- hc[which(hc$cluster == 4),]

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
    cluster_sec[i, "Mean_SD_C1"] <- paste(round(mean(cluster_1[,thisVar]), 
                                                digits = 4), "±",
                                          round(sd(cluster_1[,thisVar]), 
                                                digits = 4), sep = " ")
    cluster_sec[i, "Mean_SD_C2"] <- paste(round(mean(cluster_2[,thisVar]), 
                                                digits = 4), "±",
                                          round(sd(cluster_2[,thisVar]), 
                                                digits = 4), sep = " ")
    cluster_sec[i, "Mean_SD_C3"] <- paste(round(mean(cluster_3[,thisVar]), 
                                                digits = 4), "±",
                                          round(sd(cluster_3[,thisVar]), 
                                                digits = 4), sep = " ")
    cluster_sec[i, "p2"] <- round(t.test(cluster_1[,thisVar], 
                                         cluster_2[,thisVar])$p.value, 
                                  digits = 4)
    cluster_sec[i, "p3"] <- round(t.test(cluster_1[,thisVar], 
                                         cluster_3[,thisVar])$p.value, 
                                  digits = 4)
  }
}

# correct for multiple testing
cluster_sec$pfdr <- p.adjust(cluster_sec$p)
cluster_sec$pfdr2 <- p.adjust(cluster_sec$p2)
cluster_sec$pfdr3 <- p.adjust(cluster_sec$p3)

# ADDING HMO CLUSTERING ---------------------------------------------------------

# PREPARE THE DATA SET 
# drop unnecessary columns
rownames(meta) <- paste0("X", meta$dyad_id)

#drop meta data
meta_hmo_trim <- dplyr::select(meta, c("mom_age_at_birth", "SES", "baby_gender_cat", "mode_of_delivery_cat", "age_in_days", "pred_bf", "breastfeedings_continuous", "mom_BMI",
"Secretor", "Diversity", "Sia", "Fuc", "x2FL_nmol_ml","x3FL_nmol_ml","LNnT_nmol_ml",
"x3SL_nmol_ml","DFLac_nmol_ml","x6SL_nmol_ml","LNT_nmol_ml","LNFP_I_nmol_ml",
"LNFP_II_nmol_ml","LNFP_III_nmol_ml","LSTb_nmol_ml","LSTc_nmol_ml","DFLNT_nmol_ml",
"LNH_nmol_ml", "DSLNT_nmol_ml","FLNH_nmol_ml","DFLNH_nmol_ml","FDSLNH_nmol_ml",
"DSLNH_nmol_ml","SUM_nmol_ml", "Sia_nmol_ml","Fuc_nmol_ml","x2FL_ug_ml",
 "x3FL_ug_ml","LNnT_ug_ml","x3SL_ug_ml", "DFLac_ug_ml", "x6SL_ug_ml", 
"LNT_ug_ml", "LNFP_I_ug_ml", "breast_milk_time_hrs"))

#remove NAs from meta_hmo_trim
meta_hmo_trim <- meta_hmo_trim[complete.cases(meta_hmo_trim), ] #221 - 207 = 14 cases removed

# make a second miRNA cpm object with only secretors
table(meta_hmo_trim$Secretor)
#No Yes  No Yes 
#59 466  24 183 

secretors_hmo <- meta$dyad_id[which(meta$Secretor == "Yes")]

# only get hmo of secretors to test whether clustering is better
meta_hmo_sec <- meta[which(meta$dyad_id %in% secretors_hmo),]
meta_hmo_sec <- dplyr::select(meta_hmo_sec, c("mom_age_at_birth", "SES", "baby_gender_cat", "mode_of_delivery_cat", "age_in_days", "pred_bf", "breastfeedings_continuous", "mom_BMI",
"Secretor", "Diversity", "Sia", "Fuc", "x2FL_nmol_ml","x3FL_nmol_ml","LNnT_nmol_ml",
 "x3SL_nmol_ml","DFLac_nmol_ml","x6SL_nmol_ml","LNT_nmol_ml","LNFP_I_nmol_ml",
 "LNFP_II_nmol_ml","LNFP_III_nmol_ml","LSTb_nmol_ml","LSTc_nmol_ml","DFLNT_nmol_ml",
 "LNH_nmol_ml", "DSLNT_nmol_ml","FLNH_nmol_ml","DFLNH_nmol_ml","FDSLNH_nmol_ml",
"DSLNH_nmol_ml","SUM_nmol_ml", "Sia_nmol_ml","Fuc_nmol_ml","x2FL_ug_ml",
"x3FL_ug_ml","LNnT_ug_ml","x3SL_ug_ml", "DFLac_ug_ml", "x6SL_ug_ml", 
"LNT_ug_ml", "LNFP_I_ug_ml", "breast_milk_time_hrs"))
#183 of 221 with hmo observations were secretors

# drop the extra variables
meta_hmo_trim_ID <- meta_hmo_trim[,13:42] #207
hmo_sec_trim_ID <- meta_hmo_sec[,13:42] #183

#make numeric
meta_hmo_trim <- sapply(meta_hmo_trim_ID, as.numeric)
hmo_sec_trim <- sapply(hmo_sec_trim_ID, as.numeric)

# scale 
meta_hmo_trim_scaled <- data.frame(scale(meta_hmo_trim))
hmo_sec_trim_scaled <- data.frame(scale(hmo_sec_trim))

# now calculate the distance matrix
dist_mat <- dist(meta_hmo_trim_scaled, method = 'euclidean')
dist_mat_sec <- dist(hmo_sec_trim_scaled, method = "euclidean")
#Why Euclidian?

# cluster
hclust_avg <- hclust(dist_mat, method = "complete")
hclust_sec <- hclust(dist_mat_sec, method = "complete")
#revise methods complete 

# make a plot
dend <- as.dendrogram(hclust_avg)
dend_sec <- as.dendrogram(hclust_sec)

dend %>% plot

dend_sec %>% plot 

#change the color of the branches to correspond with the clusters

#CHANGED TO SEE WHAT THIS LOOKS LIKE CLUSTERING BY HMO
dend %>% set("branches_k_color", k = 9) %>%
  set("labels_cex", 0) %>%
  plot()

dend_sec %>% set("branches_k_color", k = 9) %>%
  set("labels_cex", 0) %>%
  plot()

plot(hclust_avg, cex = 0.4)

plot(hclust_sec, cex = 0.4)

# use the "Elbow" method to estimate the appropriate number of clusters
set.seed(78)
# Remove columns with NA values
miRNA_cpm_scaled <- miRNA_cpm_scaled[, colSums(is.na(miRNA_cpm_scaled)) == 0]
fviz_nbclust(miRNA_cpm_scaled, kmeans, method = "wss", k.max = 10)
fviz_nbclust(miRNA_sec_trim, kmeans, method = "wss", k.max = 10)
# 3 or 4 clusters,  based on this

# we should use 9 clusters based on the figure, but that results in several
# groups having only 1 person in them - try 3 and see how that looks
# with three, there is still only one individual in the third cluster (X166)
# try 2

hclust <- cutree(hclust_avg, 6)
names(hclust) <- miRNA_cpm$dyad_id
hclust <- as.data.frame(hclust)

hclust_sec2 <- cutree(hclust_sec, 4)
names(hclust) <- miRNA_sec$dyad_id

# get some info about the clustering
table(hclust)
#1   2   3   4   5   6 
#81 101  15   2   1   7 
table(hclust_sec2)
hclust_sec2
#1  2  3  4 
#99 76  1  7

# Wanting to see HMO concentration by cluster group:
cluster_assignments <- cutree(hclust_avg, k = 9)
concentrations <- meta_hmo_trim$x2FL_nmol_ml
# Create a data frame with cluster assignments and concentrations
cluster_data <- data.frame(cluster = cluster_assignments, concentrations = concentrations)

#take 2:
# Cut the dendrograms into clusters
num_clusters <- 8  # Specify the desired number of clusters
cluster_assignments_avg <- cutree(hclust_avg, k = num_clusters)
cluster_assignments_sec <- cutree(hclust_sec, k = num_clusters)

# Get the concentrations for each group
concentrations <- meta_hmo_trim_ID  # Assuming this contains the concentrations data
concentrations_sec <- hmo_sec_trim_ID  # Assuming this contains the concentrations data for secretors

# Create a data frame with cluster assignments and concentrations for the first dendrogram
cluster_data_avg <- data.frame(cluster = cluster_assignments_avg, concentrations)

# Create a data frame with cluster assignments and concentrations for the second dendrogram
cluster_data_sec <- data.frame(cluster = cluster_assignments_sec, concentrations_sec)

# Function to get the concentrations of specific variables for each cluster group
get_concentrations <- function(cluster_data, variables) {
  cluster_groups <- unique(cluster_data$cluster)
  
  result <- list()
  
  for (cluster_group in cluster_groups) {
    cluster_subset <- subset(cluster_data, cluster == cluster_group)
    concentrations_subset <- cluster_subset[, variables, drop = FALSE]
    result[[as.character(cluster_group)]] <- concentrations_subset
  }
  
  return(result)
}

# Specify the variables you want to examine
variables_of_interest <- c("cluster", "x2FL_nmol_ml", "x3FL_nmol_ml", "LNnT_nmol_ml",
                           "x3SL_nmol_ml", "DFLac_nmol_ml", "x6SL_nmol_ml", "LNT_nmol_ml",
                           "LNFP_I_nmol_ml", "LNFP_II_nmol_ml", "LNFP_III_nmol_ml",
                           "LSTb_nmol_ml", "LSTc_nmol_ml", "DFLNT_nmol_ml",
                           "LNH_nmol_ml", "DSLNT_nmol_ml", "FLNH_nmol_ml", "DFLNH_nmol_ml",
                           "FDSLNH_nmol_ml", "DSLNH_nmol_ml", "SUM_nmol_ml", "Sia_nmol_ml",
                           "Fuc_nmol_ml", "x2FL_ug_ml", "x3FL_ug_ml", "LNnT_ug_ml",
                           "x3SL_ug_ml", "DFLac_ug_ml", "x6SL_ug_ml", "LNT_ug_ml",
                           "LNFP_I_ug_ml")

# Get the concentrations for each cluster group in the only secretors dendrogram
concentrations_sec <- get_concentrations(cluster_data_sec, variables_of_interest)

# Get the concentrations for each cluster group in the first dendrogram
concentrations_avg <- get_concentrations(cluster_data_avg, variables_of_interest)

# Print the concentrations for each group in the all dendrogram
for (cluster_group in names(concentrations_avg)) {
  cat("Cluster Group", cluster_group, "in the first dendrogram:\n")
  print(concentrations_avg[[cluster_group]])
  cat("\n")
}

# Print the concentrations for each group in the only secretors dendrogram
for (cluster_group in names(concentrations_sec)) {
  cat("Cluster Group", cluster_group, "in the second dendrogram:\n")
  print(concentrations_sec[[cluster_group]])
  cat("\n")
}

# Looking at means
mean(cluster_1$x3SL_ug_ml)
mean(cluster_2$x3SL_ug_ml)
mean(cluster_3$x3SL_ug_ml)

#make data frame for the loop:
hmo_conc <- data.frame(matrix(ncol = 4))
colnames(hmo_conc) <- c("hmo", "cluster1", "cluster2", "cluster3")
#
for (hmo in variables_of_interest) {
  hmo_conc[hmo, "cluster1"] <- mean(cluster_1[,hmo], na.rm = T)
  hmo_conc[hmo, "cluster2"] <- mean(cluster_2[,hmo], na.rm = T)
  hmo_conc[hmo, "cluster3"] <- mean(cluster_3[,hmo], na.rm = T)
}
# ---------------------------------------------------------------------------
# Look at miRNA by HMO clustering ####
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

temp <- data.frame(t(summary(miRNA.pca)$importance))
temp$PC <- row.names(temp)

hc <- meta[which(!is.na(meta$cluster)),]

merged_data <- merge(meta_miRNA, hc$cluster, by = "dyad_id", all.x = TRUE)

ggbiplot(miRNA.pca, group = hc$cluster, ellipse = T, 
         var.axes = F) +
  theme_minimal() +
  labs(color = "cluster") +
  ylim(-4,2.5) + xlim(-3,3) +
  ylab("PC2 (% explained var.)") + xlab("PC1 (% explained var.)")
dev.off()

