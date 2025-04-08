import torch
import torch.nn as nn
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
import sys


# Dataset Class
class GeneExpressionDataset(Dataset):
    def __init__(self, features, labels=None):
        self.features = torch.tensor(features, dtype=torch.float32)
        self.labels = None if labels is None else torch.tensor(labels, dtype=torch.long)

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        if self.labels is not None:
            return self.features[idx], self.labels[idx]
        return self.features[idx]


# Transformer model definition
class TransformerModel(nn.Module):
    def __init__(self, input_dim, num_classes, num_heads, num_layers, d_model, dropout):
        super(TransformerModel, self).__init__()
        self.input_linear = nn.Linear(input_dim, d_model)
        self.layer_norm = nn.LayerNorm(d_model)
        self.encoder_layer = nn.TransformerEncoderLayer(d_model=d_model, nhead=num_heads, dropout=dropout)
        self.transformer_encoder = nn.TransformerEncoder(self.encoder_layer, num_layers=num_layers)
        self.classifier = nn.Linear(d_model, num_classes)

    def forward(self, x):
        x = self.input_linear(x)
        x = self.layer_norm(x)
        x = self.transformer_encoder(x)
        x = self.classifier(x)
        return x


def train_model(model, train_loader, optimizer, criterion, device, num_epochs=10):
    model.train()
    for epoch in range(num_epochs):
        for inputs, labels in train_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()


def predict_full_dataset(model, data_loader, device):
    model.eval()
    predictions = []
    with torch.no_grad():
        for inputs in data_loader:
            inputs = inputs.to(device)
            outputs = model(inputs)
            _, predicted = torch.max(outputs, 1)
            predictions.extend(predicted.cpu().numpy())
    return predictions


if __name__ == "__main__":
    # Get parameters from the command line
    args = sys.argv[1:]
    # Parameter analysis
    output_predictions_path = args[0]
    num_heads = int(args[1])
    num_layers = int(args[2])
    d_model = int(args[3])
    dropout = float(args[4])

    print(f"Output Path: {output_predictions_path}")
    print(f"Transformer parameter: num_heads={num_heads}, num_layers={num_layers}, d_model={d_model}, dropout={dropout}")

    # Data Load
    input_gene_expression_path = './scMSCF-main/preprocessing/output_data/gene_expression_afterprocess.csv'
    input_confidence_cells_path = './scMSCF-main/metaclustering/output_data/TopConfidence_cells.csv'
    gene_expression_data = pd.read_csv(input_gene_expression_path, index_col=0).transpose()
    confidence_cells = pd.read_csv(input_confidence_cells_path)
    train_cell_names = confidence_cells['cell_name'].values
    y_train = confidence_cells['final_cluster'].values - 1

    # Data alignment
    common_cells = gene_expression_data.index.intersection(train_cell_names)
    if len(common_cells) == 0:
        raise ValueError("No common cells between gene_expression_data and train_cell_names!")
    X_train = gene_expression_data.loc[common_cells].values
    train_cell_names = common_cells

    # Data preprocessing
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)

    input_dim = gene_expression_data.shape[1]
    num_classes = len(np.unique(y_train))

    # Model parameter debugging
    print(f"Model parameters: input_dim={input_dim}, num_classes={num_classes}")

    # Initialize model and training
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = TransformerModel(input_dim, num_classes, num_heads, num_layers, d_model, dropout).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.CrossEntropyLoss()

    train_dataset = GeneExpressionDataset(X_train, y_train)
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    train_model(model, train_loader, optimizer, criterion, device, num_epochs=10)

    # Full dataset prediction
    scaled_data = scaler.transform(gene_expression_data.values)
    full_dataset = GeneExpressionDataset(scaled_data)
    full_loader = DataLoader(full_dataset, batch_size=32)
    predictions = predict_full_dataset(model, full_loader, device)

    # Save forecast results
    predicted_clusters_df = pd.DataFrame({
        'Cell Name': gene_expression_data.index,
        'Predicted Cluster': predictions
    })
    predicted_clusters_df.to_csv(output_predictions_path, index=False)


