wMetaC_with_confidence <- function(nC, hmethod = hmethod, enN.cluster = enN.cluster, minN.cluster = minN.cluster, maxN.cluster = maxN.cluster, sil.thre = sil.thre, height.Ntimes, confidence_ratio = confidence_ratio) {
  # Ensure that the number of columns of nC is greater than or equal to 1
  if (ncol(nC) < 1) {
    stop("The number of nC columns must be greater than or equal to 1")
  }

  N = nrow(nC)
  C = ncol(nC)  # Get the number of columns

  # Ensure that the column index is within a reasonable range
  if (C < 1) {
    stop("The number of nC columns must be greater than or equal to 1")
  }

  # Calculate the similarity matrix of each clustering scheme
  AA = Reduce("+", apply(nC, 2, getA))
  AA = AA / C

  # Calculate the index and weight of non-zero elements
  indA = which(AA != 0, arr.ind = TRUE)
  nd = vapply(AA[indA], function(x) x * (1 - x), numeric(1))

  newAA = sparseMatrix(i = indA[, 1], j = indA[, 2], x = nd, dims = c(N, N))

  # Calculate the weight of each data point
  w0 = 4 / N * rowSums(newAA)
  e = 0.01
  w1 = (w0 + e) / (1 + e)

  # Modify Label
  x = as.vector(sapply(1:C, function(i) {
    paste(nC[, i], "_", i, sep = "")
  }))
  newnC <- matrix(x, nrow = N, byrow = FALSE)

  # Construct similarity matrix between clusters
  R = unique(x)
  allC = length(R)
  cb = combn(allC, 2)
  alls = apply(cb, 2, getss, R = R, x = x, w1 = w1)
  S0 = sparseMatrix(i = cb[1, ], j = cb[2, ], x = alls, dims = c(allC, allC))
  S = S0 + t(S0) + diag(allC)

  # Select the optimal number of clusters according to the contour coefficient
  hres = get_opt_hclust(S, hmethod, N.cluster = enN.cluster, minN.cluster, maxN.cluster, sil.thre, height.Ntimes)

  tf = hres$f
  v = hres$v
  cat("Number of clusters before meta cluster:", hres$optN.cluster, "\n")

  newnC[] <- vapply(newnC, function(q) tf[match(q, R)], numeric(1))

  finalC = apply(newnC, 1, function(d) names(sort(table(d), decreasing = TRUE)[1]))

  N.cluster = length(unique(finalC))
  perc = 0.3  # Modify voting scheme proportion parameter

  if(N.cluster == 1){
    finalC = apply(newnC, 1, function(d){
      x = sort(table(d), decreasing = TRUE)[1:2]
      n0 = length(x[1])
      if(x[2] >= n0 * perc){
        y = names(x[2])
      }else{
        y = names(x[1])
      }
      return(y)
    })
    N.cluster = length(unique(finalC))
  }

  cat("The optimal number of meta clusters:", N.cluster, "\n")

  
  file_path <- "./scMSCF-main/metaclustering/output_data/Voting_cluster.csv"
  write.csv(newnC, file = file_path, row.names = FALSE)
  cat("CSV file successfully saved to:", file_path, "\n")
  
  
  file_path <- "./scMSCF-main/metaclustering/output_data/Voting_cluster.csv"
  newnC <- read.csv(file = file_path, header = TRUE, encoding='utf-8')
  RP_times <- ncol(newnC)

  num_clusters <- ncol(newnC)
  class(num_clusters)
  num_labels <- length(unique(unlist(newnC)))

  new_columns <- paste("Cluster_", 1:enN.cluster, sep = "")
  label_columns <- colnames(newnC)
  result_df <- data.frame(matrix(0, nrow = nrow(newnC), ncol = enN.cluster))
  colnames(result_df) <- paste0("Cluster_", 1:enN.cluster)

  for (i in 1:nrow(newnC)) {
    for (cluster_num in 1:enN.cluster) {
      cluster_column <- paste0("Cluster_", cluster_num)
      count <- sum(newnC[i, label_columns] == cluster_num)
      result_df[i, cluster_column] <- count / RP_times
    }

    max_cluster <- which.max(result_df[i, ])
    newnC[i, "final_cluster"] <- max_cluster
  }

  newnC <- cbind(newnC, result_df)
  pca_data <- read.csv("./scMSCF-main/clustering/output_data/multi_clusts_combined.csv", encoding='utf-8')
  adjusted_data <- cbind(newnC, cell_name = pca_data[, 1])
  write.csv(adjusted_data, "./scMSCF-main/metaclustering/output_data/adjusted_Voting_cluster_with_cells.csv", row.names = FALSE)
  
  adjusted_file_path <- "./scMSCF-main/metaclustering/output_data/adjusted_Voting_cluster_with_cells.csv"
  adjusted_nC <- read.csv(file = adjusted_file_path, header = TRUE, encoding='utf-8')

  percent1_df <- data.frame(matrix(0, nrow = 0, ncol = ncol(adjusted_nC)))
  colnames(percent1_df) <- colnames(adjusted_nC)

  for (cluster_num in 1:enN.cluster) {
    cluster_col <- paste0("Cluster_", cluster_num)
    cluster_data <- adjusted_nC[adjusted_nC$final_cluster == cluster_num, ]

    probabilities <- cluster_data[, cluster_col] / sum(cluster_data[, cluster_col])

    top1_percent <- cluster_data[order(-probabilities), ][1:round(confidence_ratio * nrow(cluster_data)), ]

    percent1_df <- rbind(percent1_df, top1_percent)
  }

  TOP_percent_path <- "./scMSCF-main/metaclustering/output_data/TopConfidence_cells.csv"
  write.csv(percent1_df, file = TOP_percent_path, row.names = FALSE)
  cat("CSV file successfully saved to:", TOP_percent_path, "\n")

  uC = unique(finalC)
  y0 = apply(newnC, 1, function(q){
    t = rep(0, N.cluster)
    for(i in c(1:N.cluster)){
      t[i] = length(which(q %in% uC[i]))
    }
    return(t)
  })
  y0 = t(y0)

  x0 = matrix(0, nrow = N, ncol = N.cluster)
  tw = 0.5

  for(i in 1:N){
    xind = which(finalC[i]==uC)
    x0[i, xind] = 1
    allind = which(y0[i,]!=0)
    diffind = setdiff(allind, xind)
    if(length(diffind) != 0){
      x0[i, diffind] = tw * y0[i, diffind] / y0[i, xind]
    }
  }

  x0 = apply(newnC, 1, function(x) {
    t = rep(0, N.cluster)
    for (i in c(1:N.cluster)) {
      t[i] = length(which(x %in% i))
    }
    return(t)
  })
  x0 = t(x0)  # transposition

  out = list()
  out$finalC = finalC
  out$x0 = x0
  return(out)
}

