#!/usr/bin/env python

import pandas as pd
import glob
import sys
import csv
import numpy as np

accession=sys.argv[1]

# Use glob to match files with the specified wildcard
gtf_files = glob.glob(f'../../subdirs/{accession}/*_from_pipeline.gtf')

# Use glob to match files with the specified wildcard
gtf_files = glob.glob(f'../../subdirs/{accession}/*_from_pipeline.gtf')

# Define the sorting function as shown in the previous response
def gtf_sort_key(line):
    fields = line.strip().split('\t')
    seqname = fields[0]
    start = int(fields[3])
    return (seqname, start)

for gtf_file in gtf_files:
    # Create the sorted output file name
    output_file = gtf_file.replace('_from_pipeline.gtf', '_sorted.gtf')

    # Read the GTF file into a list of lines
    with open(gtf_file, 'r') as gtf_file:
        gtf_lines = gtf_file.readlines()

    # Sort the GTF lines using the custom sorting key
    sorted_gtf_lines = sorted(gtf_lines[1:], key=gtf_sort_key)

    # Write the sorted GTF lines to the output file
    with open(output_file, 'w') as sorted_file:
        # Write the header line first
        sorted_file.write(gtf_lines[0])

        # Write the sorted GTF lines
        for line in sorted_gtf_lines:
            sorted_file.write(line)

    # Create a DataFrame from the sorted GTF lines
    gtf_df = pd.read_csv(output_file, sep='\t', header=None)
# Use glob to match files with the specified wildcard
gtf_files = glob.glob(f'../../subdirs/{accession}/*_sorted.gtf')

# Initialize an empty DataFrame to store cluster information
clusters_found = pd.DataFrame(columns=['cluster_number', 'start_cl', 'stop_cl', 'chromosome', 'num_cl'])

#Determine cluster size
max_skip=2000000

# Function to check if two genes are part of the same cluster
def is_same_cluster(row1, row2):
    return row1['seqname'] == row2['seqname'] and abs(row1['start'] - row2['start']) <= max_skip

# Iterate through the GTF files
for gtf_file in gtf_files:
    # Read the GTF file with tab and semi-colon delimiters
    gtf_in = pd.read_csv(gtf_file, sep='\t', header=None, engine='python')
    gtf_in.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

    # Initialize variables to track cluster information for this file
    current_cluster = None
    start_cl = None
    stop_cl = None
    cluster_num = -1
    num_cl = 0

    gtf_in['gap'] = "NaN"
    gtf_in["cluster_id"] = "NaN"

    # Iterate through the rows of the GTF DataFrame
    for i in range(len(gtf_in)):
        row = gtf_in.iloc[i]
        cluster_temp=0

        # If no cluster is started yet, check if this row can start a new cluster
        if current_cluster is None:
            for j in range(i + 1, len(gtf_in)):
                next_row = gtf_in.iloc[j]
                gtf_in['gap'][i]="NaN"
                if is_same_cluster(row, next_row):
                    current_cluster = cluster_num
                    start_cl = row['start']
                    num_cl = 0
                    break

            # If a cluster has already started but is incomplete, check if the next gene belongs to the same cluster
        if current_cluster is not None:
            if i == len(gtf_in) - 1 or not is_same_cluster(row, gtf_in.iloc[i + 1]):
                stop_cl = row['end']
                current_cluster = None
                cluster_num += 1
                num_cl += 1
                gtf_in["cluster_id"][i]=cluster_num
#                 row1_start=gtf_in.iloc[i-1]['start']
#                 row2_start=gtf_in.iloc[i]['start']
                # Get the 'seqname' of the first entry in the cluster for the 'chromosome' column
                chromosome = gtf_in.iloc[i]['seqname']
                clusters_found = clusters_found.append({'cluster_number': cluster_num, 'start_cl': start_cl, 'stop_cl': stop_cl, 'num_cl': num_cl, 'chromosome': chromosome}, ignore_index=True)
#                 gtf_in['gap'][i] = abs(row1_start - row2_start)

            elif is_same_cluster(row, gtf_in.iloc[i + 1]):
                num_cl += 1
                row1_start=gtf_in.iloc[i]['start']
                row2_start=gtf_in.iloc[i+1]['start']
                gtf_in['gap'][i] = abs(row1_start - row2_start)
                # print(abs(row1_start - row2_start))
                gtf_in["cluster_id"][i]=cluster_num + 1

columns_to_convert = ['cluster_number','num_cl','start_cl','stop_cl']
clusters_found[columns_to_convert] = clusters_found[columns_to_convert].astype(int)

clusters_found['size'] = abs(clusters_found['stop_cl'] - clusters_found['start_cl'])
clusters_found['density'] = (clusters_found['size'] / clusters_found['num_cl']).apply(lambda x: round(x, 2))

# Initialize variables to track the cluster with the largest num_cl
largest_cluster_num = None
largest_num_cl = 0

# Iterate through the clusters_found DataFrame
for index, row in clusters_found.iterrows():
    if row['num_cl'] > largest_num_cl:
        largest_num_cl = row['num_cl']
        largest_cluster_num = row['cluster_number']


# Find the density of the largest cluster
if largest_cluster_num is not None:
    largest_cluster = clusters_found[clusters_found['cluster_number'] == largest_cluster_num]
    spacing_big = largest_cluster['density'].values[0]
else:
    spacing_big = None

# Iterate through the clusters_found DataFrame
for index, row in clusters_found.iterrows():
    start_cl_value = row['start_cl']
    stop_cl_value = row['stop_cl']


    # Locate the first row where start_cl matches in gtf_in
    j = (gtf_in.index[gtf_in['start'] == start_cl_value].min()).astype(float)

    # Locate the first row where stop_cl matches in gtf_in
    k = (gtf_in.index[gtf_in['end'] == stop_cl_value].min()).astype(float)

    median_gap = np.median(gtf_in['gap'].loc[j:k].apply(pd.to_numeric, errors='coerce').dropna())
    clusters_found.at[index, 'median'] = median_gap

#fraction of genes in clusters of 2 or more
total_genes = len(gtf_in)  # Total number of genes in the GTF file
clusters_with_2_or_more_genes = clusters_found[clusters_found['num_cl'] >= 2]
genes_in_clusters_of_2_or_more = clusters_with_2_or_more_genes['num_cl'].sum()

fraction_genes_in_clusters_of_2_or_more = genes_in_clusters_of_2_or_more / total_genes


# Calculate total_len (sum of all values of size)
total_len = clusters_found['size'].sum()

# Calculate total_num (sum of all values of num_cl)
total_num = clusters_found['num_cl'].sum()

# Calculate total_space (total_len divided by total_num)
total_space = round(total_len / total_num, 2)

# Iterate through the rows of gtf_in
for i, row in gtf_in.iterrows():
    cluster_id_value = row['cluster_id']

    # Check if gtf_in["gap"][i] is a numeric type (float or int)
    if pd.api.types.is_numeric_dtype(row['gap']):
        # Find the corresponding row in clusters_found
        cluster_row = clusters_found[clusters_found['cluster_number'] == cluster_id_value]

        if not cluster_row.empty:
            median_value = cluster_row.iloc[0]['median']

            # Check if gtf_in["gap"][i] > 5 * clusters_found["median"]
            if row['gap'] > 20 * median_value:
                gtf_in.at[i, 'break_me'] = True
            else:
                gtf_in.at[i, 'break_me'] = False
        else:
            # Handle the case where no matching cluster_number is found in clusters_found
            gtf_in.at[i, 'break_me'] = 'N/A'  # or any other appropriate value
    else:
        # Skip rows where "gap" is not a numeric type
        gtf_in.at[i, 'break_me'] = 'N/A'  # or any other appropriate value

# Create gtf_break DataFrame
gtf_break = gtf_in[gtf_in['break_me'] == True].copy()

# Create a list of unique cluster_number values from gtf_break
cluster_break = gtf_break['cluster_id'].unique().tolist()

# Create a new DataFrame gtf_out with the same columns as gtf_in
gtf_out = gtf_in.copy()

# Initialize variables
temp = 0
previous_cluster_id = None

# Iterate through the rows of gtf_in
for i in range(len(gtf_in)):
    current_row = gtf_in.iloc[i]
    break_me = gtf_in['break_me'].iloc[i-1]
    current_cluster_id = current_row['cluster_id']

    if current_cluster_id != previous_cluster_id:
        # Reset temp to 0 when cluster_id changes
        temp = 0

    if (break_me and break_me != 'N/A') or (break_me and temp == 0):
        if i < len(gtf_in) - 1 and gtf_in.iloc[i]['break_me']:
            # If both consecutive rows have break_me=True
            j = current_cluster_id
            gtf_out.at[i, 'cluster_id'] = f"{j}.NaN"
        else:
            # If this row has break_me=True but the next is False or it's the last row
            temp += 0.01
            j = current_cluster_id
            j_str=float(j)
            temp_str=float(f"{temp:.2f}")
            temp_full =j_str + temp_str
            gtf_out.at[i, 'cluster_id'] = str(f"{temp_full:.2f}")
            # print(temp_full)
            #gtf_out.at[i, 'cluster_id'] = f"{j_str + temp_str}"
    else:
        j = current_cluster_id
        gtf_out.at[i, 'cluster_id'] = f"{j + temp:.2f}"

    # Update previous_cluster_id for the next iteration
    previous_cluster_id = current_cluster_id

#Extract unique float values from the "cluster_id" column
temp_list = gtf_out['cluster_id'].apply(lambda x: f"{float(x):.2f}" if isinstance(x, (float, str)) and x.replace('.', '', 1).isdigit() else None).dropna().unique()

#Create the clusters_split DataFrame with TBD values
clusters_split = pd.DataFrame({'cluster_number': temp_list, 'start_cl': 'TBD', 'stop_cl': 'TBD', 'chromosome': 'TBD','num_cl': 0})

#Set the "cluster_number" column in clusters_split
clusters_split['cluster_number'] = temp_list

#Find the minimum and maximum values of "start" and "stop" for each cluster_number
for index, row in clusters_split.iterrows():
    cluster_number = row['cluster_number']
    mask = gtf_out['cluster_id'] == cluster_number

    if mask.any():
        min_start = gtf_out.loc[mask, 'start'].min()
        max_stop = gtf_out.loc[mask, 'start'].max()

        # Find the seqname corresponding to max_stop
        max_stop_seqname = gtf_out.loc[mask & (gtf_out['start'] == max_stop), 'seqname'].values[0]


        # Count the number of values in gtf_out with the same cluster_id
        num_cl = mask.sum()

        clusters_split.at[index, 'start_cl'] = min_start
        clusters_split.at[index, 'stop_cl'] = max_stop
        clusters_split.at[index, 'num_cl'] = num_cl

         # Set chromosome to max_stop_seqname
        clusters_split.at[index, 'chromosome'] = max_stop_seqname

clusters_split['size'] = abs(clusters_split['stop_cl'] - clusters_split['start_cl'])
clusters_split['density'] = (clusters_split['size'] / clusters_split['num_cl']).apply(lambda x: round(x, 2))


# Calculate total_len (sum of all values of size)
total_len_v2 = clusters_split['size'].sum()

# Calculate total_num (sum of all values of num_cl)
total_num_v2 = clusters_split['num_cl'].sum()

# Fraction in new clusters
fr_cl = total_num_v2 / total_genes

# Calculate total_space (total_len divided by total_num)
total_space_v2 = round(total_len_v2 / total_num_v2, 2)

# Save the clusters_found DataFrame to a CSV file
# clusters_found.to_csv('/lab/wengpj01/vertebrate_pipeline/subdirs/GCA_028390025.1/clusters_found.csv', index=False)
clusters_split.to_csv(f'../../subdirs/{accession}/clusters_found_split_2012.csv', index=False)

# Open the CSV file for writing
csv_filename = '../../results/coordinate_analysis/summary_clusters_median_1012.csv'  # Change the filename as needed
with open(csv_filename, 'a', newline='') as csv_file:
    writer = csv.writer(csv_file)

    # Write the header row if the file is empty
    if csv_file.tell() == 0:
        writer.writerow(['Accession', 'Number of Genes', 'Number of Clusters Old', 'Number of Clusters','Fraction Clustered','Total Len(kb)', 'Total Space (kb/receptor)'])

    # Write the data row
    writer.writerow([accession, total_genes, len(clusters_found), len(clusters_split), fr_cl,total_len_v2, total_space_v2])
