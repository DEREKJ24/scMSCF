import sys
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from joblib import Parallel, delayed
import multiprocessing

def process_dimension(dim, feature_matrix, optimal_k_values):
    # PCA dimension reduction
    pca = PCA(n_components=dim)
    reduced_data = pca.fit_transform(feature_matrix.T)

    results = []
    for k in optimal_k_values:
        kmeans = KMeans(n_clusters=k, random_state=40)
        cluster_labels = kmeans.fit_predict(reduced_data)
        column_name = f'cluster_dim{dim}_k{k}'
        results.append((column_name, cluster_labels))

    return results

def main():
    # Extract command line parameters
    optimal_k_values_list = eval(sys.argv[1])  # Combination of K values entered by the user

    # Path Configuration
    csv_file_path = './scMSCF-main/preprocessing/output_data/gene_expression_2000genes.csv'
    result_csv_path = './scMSCF-main/clustering/output_data/multi_clusts_combined.csv'

    # Data reading
    feature_matrix = pd.read_csv(csv_file_path, index_col=0)

    # Create an empty DataFrame to store all clustering results
    combined_clusters_df = pd.DataFrame(index=feature_matrix.columns)
    combined_clusters_df.index.name = 'Barcode'

    # Process each dimension and K value in parallel
    results = Parallel(n_jobs=multiprocessing.cpu_count())(
        delayed(process_dimension)(dim, feature_matrix, optimal_k_values)
        for dim, optimal_k_values in optimal_k_values_list
    )

    # Integration results
    for result in results:
        for column_name, cluster_labels in result:
            combined_clusters_df[column_name] = cluster_labels

    # Save results to CSV
    combined_clusters_df.to_csv(result_csv_path)
    print(f"All clustering completed and labels saved to {result_csv_path}")

if __name__ == '__main__':
    main()
