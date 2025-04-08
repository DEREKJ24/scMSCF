scMCTL is an innovative clustering algorithm that seamlessly integrates multi-dimensional PCA for dimensionality reduction, K-means clustering, and a weighted ensemble meta-clustering approach. Enhanced by a self-attention-powered Transformer model, scMCTL is designed to elevate clustering performance in single-cell RNA sequencing analysis.

scMCTL employs a multi-layer dimensionality reduction strategy to construct an initial clustering framework, establishing a robust consensus clustering structure. Through a voting mechanism embedded in the meta-clustering process, high-confidence cells are selected from the initial clustering results, providing precise training labels for the Transformer model. This enables the model to capture complex dependencies in gene expression data, significantly enhancing clustering accuracy and reliability.

**Step 1: Data Preprocessing**  
  Input: Raw gene expression matrix (CSV format).  
  Output: Processed gene expression matrix and a subset of 2000 highly variable genes.  
  Script: preprocessing/preprocessing.R  

**Step 2: PCA and K-means Clustering**  
  Input: Processed gene expression data.  
  Output: Candidate clusters generated through PCA and K-means clustering.  
  Script: clustering/PCA_multiK_cluster.py  

**Step 3: Weighted Metaclustering**    
  Input: Multiple clustering results from Step 2.  
  Output: High-confidence clusters generated through weighted metaclustering.  
  Script: metaclustering/Main_wMetaC.R  

**Step 4: Transformer-Based Modeling**    
  Input: Gene expression data and high-confidence clusters.  
  Output: Final predicted clusters.  
  Script: transformer/transformer4.py  

  
![Figure 1](https://github.com/user-attachments/assets/9d4abeca-0603-4264-9b60-daffd35405f1)
