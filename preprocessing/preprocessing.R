# Get file path from command line parameters
args <- commandArgs(trailingOnly = TRUE)
raw_data_path <- args[1]
output_processed_path <- args[2]
output_variable_genes_path <- args[3]

# Load gene expression matrix
gene_expression <- read.csv(raw_data_path, row.names = 1)
gene_expression_matrix <- as.matrix(gene_expression)

# Create Seurat objects and preprocess them
library(Seurat)
seurat_object <- CreateSeuratObject(counts = gene_expression_matrix)
seurat_object <- SCTransform(seurat_object)

# Save the complete matrix processed by SCTransform
write.csv(seurat_object@assays$SCT@data, file = output_processed_path)

# Extraction and preservation of 2000 hypervariable genes
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
variable_genes <- VariableFeatures(seurat_object)
filtered_expression_matrix <- seurat_object@assays$SCT@data[variable_genes, ]
write.csv(filtered_expression_matrix, file = output_variable_genes_path)
