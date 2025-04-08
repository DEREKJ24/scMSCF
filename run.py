import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


def run_r_script(script_path, args=None):
    """General function for running R script"""
    rscript_path = r'.\Rscript.exe'  # Replace with actual path
    cmd = [rscript_path, script_path]
    if args:
        cmd += args
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8')
    if result.returncode != 0:
        print(f"Error running {script_path}: {result.stderr}")
    else:
        print(f"Successfully ran {script_path}")


def run_python_script(script_path, args=None):
    """General functions for running Python scripts"""
    cmd = ['python', script_path]
    if args:
        cmd += args
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8')
    if result.returncode != 0:
        print(f"Error running {script_path}: {result.stderr}")
    else:
        print(f"Successfully ran {script_path}")


def plot_elbow_method(data, max_k):
    """Generate elbow normal graph to help users select the best K value"""
    sse = []
    k_values = range(2, max_k + 1)
    for k in k_values:
        kmeans = KMeans(n_clusters=k, random_state=40)
        kmeans.fit(data)
        sse.append(kmeans.inertia_)
    plt.figure(figsize=(10, 6))
    plt.plot(k_values, sse, marker='o')
    plt.xlabel('Number of clusters (K)')
    plt.ylabel('Sum of Squared Errors (SSE)')
    plt.title('Elbow Method for Determining Optimal K')
    plt.show()
    optimal_k_values = input(
        "Please input two optimal K values as observed from the elbow plot, separated by a comma: ").split(',')
    return [int(k) for k in optimal_k_values]


def step_1_preprocessing(raw_data_path, output_processed_path, output_variable_genes_path, preprocessing_script):
    """Step 1: Data processing"""
    print("Step 1: Data preprocessing start...")
    os.makedirs(os.path.dirname(output_processed_path), exist_ok=True)
    os.makedirs(os.path.dirname(output_variable_genes_path), exist_ok=True)
    run_r_script(preprocessing_script, args=[raw_data_path, output_processed_path, output_variable_genes_path])
    print("Step 1: Data preprocessing completed！")


def step_2_pca_kmeans(pca_file_path, dimensions_list, max_k, clustering_script):
    """Step 2: PCA dimension reduction and K-means clustering"""
    print("Step 2: PCA dimension reduction and K-means clustering start...")
    feature_matrix = pd.read_csv(pca_file_path, index_col=0)
    optimal_k_values_list = []
    for dim in dimensions_list:
        print(f"Generating elbow plot for dimension: {dim}")
        pca = PCA(n_components=dim)
        reduced_data = pca.fit_transform(feature_matrix.T)
        optimal_k_values = plot_elbow_method(reduced_data, max_k)
        print(f"Optimal K values for dimension {dim}: {optimal_k_values}")
        optimal_k_values_list.append((dim, optimal_k_values))
    os.makedirs('./scMSCF-main/scMSCF-main/clustering/output_data', exist_ok=True)
    run_python_script(clustering_script, args=[str(optimal_k_values_list)])
    print("Step 2: PCA dimension reduction and K-means clustering are completed!")


def step_3_metaclustering(input_clusters_path, target_cluster, confidence_ratio, metaclustering_script):
    """Step 3: Weighted Metacluster"""
    print("Step 3: weighted meta cluster start...")
    run_r_script(metaclustering_script, args=[input_clusters_path, str(target_cluster), str(confidence_ratio)])
    print("Step 3: Weighted meta clustering is complete!")


def step_4_transformer(output_predictions_path, num_heads, num_layers, d_model, dropout, transformer_script):
    """Step 4: Transformer model training and prediction"""
    print("Step 4: Transformer model training and prediction start...")
    os.makedirs(os.path.dirname(output_predictions_path), exist_ok=True)
    run_python_script(transformer_script, args=[
        output_predictions_path,
        str(num_heads),
        str(num_layers),
        str(d_model),
        str(dropout)
    ])
    print("Step 4: Transformer model training and prediction completed!")


if __name__ == "__main__":
    # User configuration parameters
    raw_data_path = input("Please enter the path of the original gene expression matrix (for example: E:/scMSCF-main/preprocessing/input_data/SimCell.csv):")
    confidence_ratio = float(input("Please enter the confidence ratio according to the dataset size (range 0-1, default: 0.5): ") or 0.5)
    num_heads = int(input("Please enter the number of multi head attention heads of the Transformer model according to the dataset size (default: 8): ") or 8)
    num_layers = int(input("Please enter the number of Transformer encoder layers according to the data set size (default: 4): ") or 4)
    d_model = int(input("Please enter the d_model value of Transformer according to the data set size (default: 768): ") or 768)
    dropout = float(input("Please refer to Dropout probability (default: 0.1) based on dataset size:") or 0.1)

    # Step 1
    preprocessing_script = './scMSCF-main/preprocessing/preprocessing.R'
    output_processed_path = './scMSCF-main/preprocessing/output_data/gene_expression_afterprocess.csv'
    output_variable_genes_path = './scMSCF-main/preprocessing/output_data/gene_expression_2000genes.csv'

    # Step 2
    pca_file_path = output_variable_genes_path
    dimensions_list = [30, 35, 40, 45, 50]
    max_k = 20
    clustering_script = './scMSCF-main/clustering/PCA_multiK_cluster.py'
    # D:/GOOGLE_Download/scMSCF-main/preprocessing/input_data/SimCell.csv
    # Step 3
    input_clusters_path = './scMSCF-main/clustering/output_data/multi_clusts_combined.csv'
    metaclustering_script = './scMSCF-main/metaclustering/Main_wMetaC.R'
    target_cluster = int(input("Please enter the number of target clusters (positive integer, default: 8): ") or 8)

    # Step 4
    output_predictions_path = './scMSCF-main/transformer/output_data/predicted_clusters_full.csv'
    transformer_script = './scMSCF-main/transformer/transformer.py'

    # Execution
    step_1_preprocessing(raw_data_path, output_processed_path, output_variable_genes_path, preprocessing_script)
    step_2_pca_kmeans(pca_file_path, dimensions_list, max_k, clustering_script)
    step_3_metaclustering(input_clusters_path, target_cluster, confidence_ratio, metaclustering_script)
    step_4_transformer(output_predictions_path, num_heads, num_layers, d_model, dropout, transformer_script)
