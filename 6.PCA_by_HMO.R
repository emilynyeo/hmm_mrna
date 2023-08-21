# HEADER -----------------------------------------------------------------------
#
# TITLE:   6-PCA_BY_HMO_CLUSTERS.R
#
# PURPOSE: do a PCA analysis of miRNAs by HMO cluster
#
# DATE:    March 14, 2022
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
               readxl, factoextra, ggbiplot) 

# REad in meta data
meta <- read.csv("input/meta_clean_EY.csv")

# Do miRNA PCA: ####
miRNA_cpm <- read.csv("input/miRNA_counts_EY.csv") # reading in for colnames
# Read in PCA data from script 2
miRNA_cpm_noQC <- read_csv("input/miRNA_cpm_noQC_EY.csv")

# calculate principle components
# miRNA.pca <- prcomp(miRNA_cpm_noQC, center = T, scale = T) # non numeric error
miRNA.pca <- prcomp(miRNA_cpm_noQC[,2:210], center = T, scale = T)

# check whether we can use eigenvalue > 1 as a metric
I((summary(miRNA.pca)$sdev)^2) #this might not be NB. Not great. Too many
# eigenvalue = standard deviation^2. Too many to use those > 1.
plot(miRNA.pca)
#save the PCA results in a data frame
PCA <- data.frame(miRNA.pca$x)

# Cluster HMOs: ####

# Read in HMO clustering data from script 3
meta_hmo_trim_ID <- read_csv("input/meta_hmo_trim_EY.csv")
# Only Secretos 
hmo_sec_trim_ID <- read_csv("input/hmo_sec_trim_EY.csv")

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

hclust <- cutree(hclust_avg, 6)
names(hclust) <- miRNA_cpm$dyad_id
#hclust <- as.data.frame(hclust)

secretors <- meta$dyad_id[which(meta$Secretor == "Yes")]
miRNA_sec <- miRNA_cpm[which(miRNA_cpm$X %in% secretors),]
hclust_sec2 <- cutree(hclust_sec, 6)
names(hclust_sec2) <- miRNA_sec$dyad_id
#hclust_sec2 <- as.data.frame(hclust_sec2)

# get some info about the clustering
table(hclust)
#1   2   3   4   5   6 
#81 101  15   2   1   7 
table(hclust_sec2)
#1  2  3  4 
#99 76  1  7

#make a list of individuals in cluster 1 - cluster 3 for all: 
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

#only keep individuals in meta who have miRNA data
hc_pca <- hc[hc$dyad_id %in% row.names(PCA),] #51 in cluster 1-3 with miRNA
PCA$dyad_id <- row.names(PCA)

# visualizations for Secretor = Yes vs. No
ggbiplot(miRNA.pca, group = hc_pca$Secretor, ellipse = T, 
         var.axes = F) +
  theme_minimal() +
  labs(color = "Breast\nFeeding\nCategory") +
  ylim(-4,2.5) + xlim(-3,3) +
  ylab("PC2 (15.6% explained var.)") + xlab("PC1 (18.6% explained var.)")

#
