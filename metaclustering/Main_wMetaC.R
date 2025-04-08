# #Set CRAN image (for example, Tsinghua University image)
# options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(SHARP)
library(Matrix)
library(cluster)
library(ggplot2)


# Get file path from command line parameters
args <- commandArgs(trailingOnly = TRUE)
multi_clusts <- args[1]
clust_num <- as.integer(args[2])
confidence_ratio <- as.integer(args[3])

# Read cluster information
data <- read.csv(multi_clusts, encoding='utf-8')

# Get the column with "cluster_dim" in the column name
cluster_columns <- colnames(data)[grep("cluster_dim", colnames(data))]

# Extract these columns
nC <- data[, cluster_columns]

# Load the script containing the custom function "wMetaC_with_confidence"
source("./scMSCF-main/metaclustering/my_wMetaC.R")


# Run the wMetaC_with_confidence function
result <- wMetaC_with_confidence(nC, hmethod = "ward.D", enN.cluster = clust_num, minN.cluster = 1, maxN.cluster = 20, sil.thre = 0.5 , confidence_ratio = confidence_ratio)



