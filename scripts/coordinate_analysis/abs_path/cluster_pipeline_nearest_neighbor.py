#!/usr/bin/env python
import pandas as pd
import numpy as np
import os
import sys
import csv
import glob


#SetWD

# Specify the path to the directory you want to change to
accession = sys.argv[1]
new_directory = os.path.join("/lab/wengpj01/vertebrate_pipeline/subdirs/", accession)

print(new_directory)

# Use os.chdir() to change the working directory
os.chdir(new_directory)


# Read in the GTF file as a dataframe
coordinates = pd.read_csv(sys.argv[2], sep="\t", header=None)
# coordinates = pd.read_csv("/lab/wengpj01/vertebrate_pipeline/clusters/wood_frog.gtf", sep="\t", header=None)
coordinates.columns = ["chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

# Read in the CSV file as a dataframe with the first row as a header
clusters = pd.read_csv(sys.argv[3])

clusters = clusters.rename(columns={'cluster number': 'cluster_number'})
clusters_found = clusters
# clusters = pd.read_csv("/lab/wengpj01/vertebrate_pipeline/clusters/wood_frog_clusters.csv")

# Read in the distance matrices as pandas DataFrames
evo = pd.read_csv(sys.argv[4], sep='\t')

max_skip = int(sys.argv[5])
max_skip_str = str(max_skip)


# evo = pd.read_csv('/lab/wengpj01/vertebrate_pipeline/trees/mini_wood_frog.dist', sep='\t')
# output_phys_evo = "/lab/wengpj01/vertebrate_pipeline/trees/wood_combo.tsv"

print("Rows of coordinates:", len(coordinates))
print("Rows of evo:", len(evo))


# Get the number of rows and columns in evo
n_rows, n_cols = evo.shape

# Fill in the lower left triangle of evo
for i in range(n_rows):
    for j in range(i):
        evo.iloc[i, j] = evo.iloc[j, i]


# Initialize the "cluster_match" column in the coordinates dataframe with "TBD"
coordinates["cluster_match"] = "NaN"

# Iterate through each row in the clusters dataframe and update the "cluster_match" column
for index, cluster_row in clusters.iterrows():
    chromosome = cluster_row["chromosome"]
    start_cl = cluster_row["start_cl"]
    stop_cl = cluster_row["stop_cl"]
    cluster_number = cluster_row["cluster_number"]

    # Check if there is a match in the coordinates dataframe
    mask = (coordinates["chromosome"] == chromosome) & \
           (coordinates["start"] >= (start_cl - 10)) & \
           (coordinates["end"] <= (stop_cl + 10))

    # Update the "cluster_match" column for matching rows
    coordinates.loc[mask, "cluster_match"] = cluster_number

# Create a new DataFrame 'nearest' with the specified columns
nearest = pd.DataFrame(columns=['query', 'match', 'evo_m', "cluster_m"], index=range(evo.shape[1]))
nearest['query'] = nearest.index
nearest['match'] = "TBD"
nearest['cluster_m'] = "TBD"
nearest['evo_m'] = "TBD"

same_cluster=0
different_cluster=0
total=0

# Find the minimum non-zero values in the 'evo' DataFrame and populate 'nearest'
for i in range(evo.shape[0]):
    total += 1
    row_values = evo.iloc[i, :].values
    # Check if the row contains at least one non-zero and non-NaN value
    valid_values = row_values[(row_values != 0) & ~np.isnan(row_values)]
    if len(valid_values) == 0:
        continue  # Skip rows with all zero, all NaN, or a mix of zero and NaN values
    min_value = np.nanmin(valid_values)
    j = np.where(row_values == min_value)[0][0]
    nearest.at[i, 'match'] = j
    nearest.at[i, 'evo_m'] = min_value

    # Copy the value from 'phys' DataFrame, including NaN values
    cluster_1 = coordinates['cluster_match'][i]
    cluster_2 = coordinates['cluster_match'][j]


    if cluster_1 == cluster_2:
        nearest.at[i, 'cluster_m'] = cluster_1
        same_cluster +=1
    else:
        nearest.at[i, 'cluster_m'] = "Different"  # Ensure NaN is copied over
        different_cluster += 1


Fraction_same = same_cluster / total
print(f"Same cluster: {Fraction_same}")
# Save the 'nearest' DataFrame as a TSV file
# nearest.to_csv(output_phys_evo, sep='\t', index=False)

#fraction of genes in clusters of 2 or more
# total_genes = len(input)  # Total number of genes in the GTF file
clusters_with_2_or_more_genes = clusters_found[clusters_found['num_cl'] >= 2]
genes_in_clusters_of_2_or_more = clusters_with_2_or_more_genes['num_cl'].sum()

fraction_genes_in_clusters_of_2_or_more = round(genes_in_clusters_of_2_or_more / total,2)

# Initialize variables to track the cluster with the largest num_cl
largest_cluster_num = None
largest_num_cl = 0

# Iterate through the clusters_found DataFrame
for index, row in clusters_found.iterrows():
    if row['num_cl'] > largest_num_cl:
        largest_num_cl = row['num_cl']
        largest_cluster_num = row['cluster_number']

# # Find the spacing of the largest cluster
# if largest_cluster_num is not None:
#     largest_cluster = clusters_found[clusters_found['cluster_number'] == largest_cluster_num]
#     genes_in_largest_cluster = int(largest_cluster['num_cl'])
#     spacing_big = round(largest_cluster['spacing'].values[0]/1000,2)
# else:
#     spacing_big = None


# Calculate total_len (sum of all values of size)
total_len = round(clusters_found['size'].sum()/1000,2)

# Calculate total_num (sum of all values of num_cl)
total_num = clusters_found['num_cl'].sum()

# Calculate total_space (total_len divided by total_num)
total_space = round(total_len / total_num, 2)

# # Calculate the fraction of genes in the largest cluster
# fraction_genes_in_largest_cluster = round(genes_in_largest_cluster / total, 2)



# Open the CSV file for writing
csv_filename = glob.glob(f'/lab/wengpj01/vertebrate_pipeline/20231006_run/summary_clusters_withnearest_{max_skip_str}.csv')
print(csv_filename)
with open(csv_filename[0], 'a', newline='') as csv_file:
    writer = csv.writer(csv_file)

    # Write the header row if the file is empty
    if csv_file.tell() == 0:
        writer.writerow(['Accession', 'Number of Genes', 'Number of Clusters', 'Fraction Clustered', 'Total Len(kb)', 'Total Space (kb/receptor)', 'Nearest_in_cluster'])

    # Write the data row
    writer.writerow([accession, total, len(clusters_found), fraction_genes_in_clusters_of_2_or_more, total_len, total_space,Fraction_same])
