#!/usr/bin/env python
#SBATCH --job-name=networks # Job name
#SBATCH --mem=200gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/networks_%j.log # Standard output and error log
#SBATCH --error=logs/networks_%j.err
import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import rankdata

# Read the matrix from the CSV file
file_path = '../../results/phylogenetics/hg_p_value_v7_q2.csv'

# Assuming row 1 contains column names and column 1 contains row names
data = pd.read_csv(file_path, index_col=0)

# Convert the data to a numpy array for matrix operations
matrix = data.to_numpy()

# def benjamini_hochberg_correction(p_values):
#     # Sort p-values and get corresponding ranks
#     sorted_indices = np.argsort(p_values)
#     sorted_p_values = p_values[sorted_indices]
#     ranks = rankdata(sorted_p_values)
#
#     # Calculate Benjamini-Hochberg critical values
#     n = len(p_values)
#     bh_values = sorted_p_values * n / ranks
#
#     # Adjusted p-values
#     adjusted_p_values = np.minimum(bh_values, 1.0)
#
#     # Rearrange adjusted p-values according to original order
#     adjusted_p_values_original_order = np.empty_like(p_values)
#     adjusted_p_values_original_order[sorted_indices] = adjusted_p_values
#
#     return adjusted_p_values_original_order
# #
# # adjusted_p_values = benjamini_hochberg_correction(matrix)
# corrected_p_values = benjamini_hochberg_correction(matrix)

# Number of comparisons (matrix elements)
n_comparisons = np.prod(matrix.shape)

# Calculate Bonferroni-corrected p-values
corrected_p_values = np.minimum(matrix * n_comparisons, 1.0)
#corrected_p_values = matrix

# matrix = np.array([
#     [0, 0.1, 0.001, 0.25],
#     [0.1, 0, 0.2, 0.001],
#     [0.001, 0.2, 0, 0.5],
#     [0.5, 0.04, 0.2, 0]
# ])

# Define a threshold
threshold = 0.05

# Create a graph
G = nx.Graph()

# Add edges to the graph based on the matrix
rows, cols = np.where(corrected_p_values < threshold)
edges = zip(rows, cols)

for edge in edges:
    G.add_edge(edge[0], edge[1])

# Find connected components in the graph
clusters = list(nx.connected_components(G))

# # Create a graph
# G = nx.Graph()
#
# # Add edges to the graph based on the matrix with at least two connections
# for i in range(len(matrix)):
#     for j in range(i + 1, len(matrix)):
#         if np.sum(matrix[i] < threshold) >= 2 and np.sum(matrix[j] < threshold) >= 2:
#             if np.sum((matrix[i] < threshold) & (matrix[j] < threshold)) >= 2:
#                 G.add_edge(i, j)
#
# # Find connected components in the graph
# clusters = list(nx.connected_components(G))

# print(clusters)

output_file = '../../results/phylogenetics/clusters_output_og_v7_q2_withBH.txt'

with open(output_file, 'w') as file:
    for idx, cluster in enumerate(clusters, 1):
        file.write(f'Cluster {idx}: {cluster}\n')

# Map indices to their corresponding row and column names
row_names = data.index.tolist()
col_names = data.columns.tolist()

named_clusters = []
for cluster in clusters:
    named_cluster = {
        'rows': [row_names[i] for i in cluster],
        'columns': [col_names[i] for i in cluster]
    }
    named_clusters.append(named_cluster)

# Write named clusters to a text file
output_file = '../../results/phylogenetics/clusters_output_with_names_v7_q2_withBH.txt'

with open(output_file, 'w') as file:
    for idx, cluster in enumerate(named_clusters, 1):
        file.write(f'Cluster {idx}:\n')
        file.write(f'Rows: {cluster["rows"]}\n')
        file.write(f'Columns: {cluster["columns"]}\n\n')



# Display clusters
for idx, cluster in enumerate(clusters, 1):
    print(f'Cluster {idx}: {cluster}')
