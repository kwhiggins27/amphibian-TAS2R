#!/usr/bin/env python
#SBATCH --job-name=getqk_2fast # Job name
#SBATCH --mail-type=FAIL,END    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute    # Where to send mail
#SBATCH --mem=20gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/getqk_2fast_%j.log # Standard output and error log
#SBATCH --error=logs/getqk_2fast_%j.err
import pandas as pd
import os
accepted_range = pd.read_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/maximal_range_all.tsv', sep='\t')  #This df is notable for columns chromosome, start, stop, accession
accepted_range['id'] = (
    accepted_range['accession'] + "_" +
    accepted_range['chromosome'].astype(str) + ":" +
    accepted_range['start'].astype(str) + "_" +
    accepted_range['stop'].astype(str)
)
calc_N = pd.DataFrame(index=accepted_range["id"], columns=accepted_range["id"]).reset_index(drop=True)
calc_m = pd.DataFrame(index=accepted_range["id"], columns=accepted_range["id"]).reset_index(drop=True)
calc_n = pd.DataFrame(index=accepted_range["id"], columns=accepted_range["id"]).reset_index(drop=True)
calc_q = pd.DataFrame(index=accepted_range["id"], columns=accepted_range["id"]).reset_index(drop=True)
calc_k = pd.DataFrame(index=accepted_range["id"], columns=accepted_range["id"]).reset_index(drop=True)
columns = ['accession', 'Busco id', 'Status', 'Sequence', 'Gene Start', 'Gene End', 'Strand', 'Score', 'Length', 'OrthoDB url', 'Description']

# Iterating through pairs of accepted_range_A and accepted_range_B
for index_A, row_A in accepted_range.iterrows():
    chromosome_A = str(row_A['chromosome'])
    start_A = row_A['start']
    stop_A = row_A['stop']
    accession_A = row_A['accession'].split(':')[0]  # Extracting accession for species A

    file_path_A = f'/lab/wengpj01/vertebrate_pipeline/busco/{accession_A}/busco_output/run_vertebrata_odb10/full_table.tsv'

    if os.path.exists(file_path_A):
        try:
            species_A = pd.read_csv(file_path_A, sep='\t', skiprows=3)
            print(f"Read successful for {file_path_A}")
        except pd.errors.ParserError as e:
            print(f"Error reading {file_path_A}: {e}")
            continue

        if species_A.shape[1] != 10:
            calc_q.iloc[index_A, :] = 0
            calc_k.iloc[index_A, :] = 0
            continue
    else:
        calc_q.iloc[index_A, :] = 0
        calc_k.iloc[index_A, :] = 0
        continue
    accepted_range2 = accepted_range[index_A:]
    for index_B, row_B in accepted_range2.iterrows():
        chromosome_B = str(row_B['chromosome'])
        start_B = row_B['start']
        stop_B = row_B['stop']
        accession_B = row_B['accession'].split(':')[0]  # Extracting accession for species B

        file_path_B = f'/lab/wengpj01/vertebrate_pipeline/busco/{accession_B}/busco_output/run_vertebrata_odb10/full_table.tsv'

        # Check if files exist
        if os.path.exists(file_path_B):
            # Read species A and species B files, skipping the first 3 rows
            try:
                species_B = pd.read_csv(file_path_B, sep='\t', skiprows=3)
                print(f"Read successful for {file_path_B}")
            except pd.errors.ParserError as e:
                print(f"Error reading {file_path_B}: {e}")
                continue
            # Check if busco_output has 10 columns
            if species_A.shape[1] != 10 or species_B.shape[1] != 10:
                continue

            species_A.columns = columns[1:]
            species_B.columns = columns[1:]

            common_genes_all = set(species_A.iloc[:, 0]).intersection(species_B.iloc[:, 0])
            calc_N.iloc[index_A, index_B] = len(common_genes_all)

            #Find subset of species_A that matches this specific clusterA
             # Drop rows where 'Gene Start' is NaN

            species_A_subset = species_A[species_A['Busco id'].isin(common_genes_all)].reset_index(drop=True)  # Reset index for species_A
            species_A_subset['Sequence'] = species_A_subset['Sequence'].str.split(':').str[0]  # Reset index for accepted_range
            species_A_subset = species_A_subset.dropna(subset=['Gene Start'])
            species_A_subset['Gene Start'] = species_A_subset['Gene Start'].astype(str).str.rstrip('.0').astype(int)
            species_A_subset = species_A_subset.dropna(subset=['Gene End'])
            species_A_subset['Gene End'] = species_A_subset['Gene End'].astype(str).str.rstrip('.0').astype(int)

            # Check conditions for each row in busco_output against the current row in subset
            matchesA = species_A_subset[
                (species_A_subset['Sequence'] == chromosome_A) &
                (species_A_subset['Gene Start'] >= start_A) &
                (species_A_subset['Gene End'] >= start_A) &
                (species_A_subset['Gene Start'] <= stop_A) &
                (species_A_subset['Gene End'] <= stop_A)
                ]
            calc_k.iloc[index_A, index_B] = len(matchesA)

            #Find subset of species_B that matches specific clusterB:
            species_B_subset = species_B[species_B['Busco id'].isin(common_genes_all)].reset_index(drop=True)  # Reset index for species_A
            accepted_range_subset_B = accepted_range[accepted_range['accession']==accession_B].reset_index(drop=True)
            species_B_subset['Sequence'] = species_B_subset['Sequence'].str.split(':').str[0]  # Reset index for accepted_range

             # Drop rows where 'Gene Start' is NaN
            species_B_subset = species_B_subset.dropna(subset=['Gene Start'])
            species_B_subset['Gene Start'] = species_B_subset['Gene Start'].astype(str).str.rstrip('.0').astype(int)
            species_B_subset = species_B_subset.dropna(subset=['Gene End'])
            species_B_subset['Gene End'] = species_B_subset['Gene End'].astype(str).str.rstrip('.0').astype(int)

                # Check conditions for each row in busco_output against the current row in subset
            matchesB = species_B_subset[
                (species_B_subset['Sequence'] == chromosome_B) &
                (species_B_subset['Gene Start'] >= start_B) &
                (species_B_subset['Gene End'] >= start_B) &
                (species_B_subset['Gene Start'] <= stop_B) &
                (species_B_subset['Gene End'] <= stop_B)
                ]
            # Extracting unique Busco ids from matchesA and matchesB

            calc_m.iloc[index_A, index_B] = len(matchesB)
            calc_n.iloc[index_A, index_B] = len(common_genes_all) - len(matchesB)

            busco_ids_matchesA = set(matchesA['Busco id'])
            busco_ids_matchesB = set(matchesB['Busco id'])

            # Finding common unique Busco ids between matchesA and matchesB
            common_busco_ids = busco_ids_matchesA.intersection(busco_ids_matchesB)

            # Number of unique Busco ids common in both matchesA and matchesB
            calc_q.iloc[index_A, index_B] = len(common_busco_ids)
        else:
            calc_m.iloc[index_A, index_B] = 0  # If files don't exist, set the value to 0
            calc_n.iloc[index_A, index_B] = 0  # If files don't exist, set the value to 0
            calc_k.iloc[index_A, index_B] = 0  # If files don't exist, set the value to 0
            calc_q.iloc[index_A, index_B] = 0  # If files don't exist, set the value to 0
calc_m.to_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/calc_m7.csv')
calc_n.to_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/calc_n7.csv')
calc_k.to_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/calc_k7.csv')
calc_q.to_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/calc_q7.csv')
