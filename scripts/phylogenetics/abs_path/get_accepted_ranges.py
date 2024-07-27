#!/usr/bin/env python
import pandas as pd
import os

# Read accession numbers from the file
accession_file = '/lab/wengpj01/vertebrate_pipeline/accessions_to_keep.txt'
with open(accession_file, 'r') as f:
    accession_list = [line.strip() for line in f.readlines()]
# Initialize bed0 DataFrame
bed3 = pd.DataFrame(columns=['chromosome', 'start', 'stop', 'accession'])

for accession in accession_list:
    bed1_file = f'/lab/wengpj01/vertebrate_pipeline/subdirs/{accession}/clusters_margins_1M_v3.bed'

    # Check if the bed1 file exists
    if not os.path.exists(bed1_file):
        bed1=pd.DataFrame(columns=['chromosome', 'start', 'stop'])

    try:
        # Read bed1 as DataFrame using whitespace as delimiter
        column_names = ['chromosome', 'start', 'stop']
        bed1 = pd.read_csv(bed1_file, delim_whitespace=True, header=None, names=column_names)
        if len(bed1) == 0:
            bed1=pd.DataFrame(columns=['chromosome', 'start', 'stop'])
        # Check if the number of rows in bed1 is even
        if len(bed1) % 2 != 0:
            print(f"Number of rows in bed1 for accession {accession} is not even. Exiting...")
    except:
        bed1=pd.DataFrame(columns=['chromosome', 'start', 'stop'])
    bed0_file = f'/lab/wengpj01/vertebrate_pipeline/subdirs/{accession}/singletons_margins_1M_v3.bed'

    # # Check if the bed1 file exists
    # if not os.path.exists(bed0_file):
    #     print(f"File {bed1_file} does not exist for accession {accession}. Skipping...")
    #     continue

    try:
        # Read bed1 as DataFrame using whitespace as delimiter
        column_names = ['chromosome', 'start', 'stop']  # Specify column names
        bed0 = pd.read_csv(bed0_file, delim_whitespace=True, header=None, names=column_names)
        if len(bed0) == 0:
            bed0=pd.DataFrame(columns=['chromosome', 'start', 'stop'])
        # Check if the number of rows in bed1 is even
        if len(bed0) % 2 != 0:
            print(f"Number of rows in bed1 for accession {accession} is not even. Exiting...")
            break
    except:
        bed0=pd.DataFrame(columns=['chromosome', 'start', 'stop'])

    bed01 = pd.concat([bed1, bed0], ignore_index=True, axis=0)


    if len(bed01)==0:
        print(f"Neither singletons nor clusters found for {accession}")
        continue

    if bed01 is not None and not bed01.empty:

        # Initialize bed2 DataFrame
        bed2 = pd.DataFrame(columns=['chromosome', 'start', 'stop', 'accession'])
        for i in range(0, len(bed01), 2):
            # Calculate new rows for bed2
            new_row = {
                'chromosome': bed01.iloc[i, 0],
                'start': max(0,bed01.iloc[i, 2] - 1000000),
                'stop': bed01.iloc[i + 1, 1] + 1000000,
                'accession': accession
            }
            bed2 = bed2.append(new_row, ignore_index=True)
        # Append bed2 to bed0
        bed3 = bed3.append(bed2, ignore_index=True)
    else:
        print("Neither singletons nor clusters found")
bed4 = bed3[(bed3['chromosome'] != 0) & (bed3['start'].notnull()) & (bed3['stop'].notnull())].copy()

# Save bed0 as maximal_range.tsv
bed4.to_csv('/lab/wengpj01/vertebrate_pipeline/search_neighborhood/maximal_range_all.tsv', sep='\t', index=False)

print(bed4[bed4["accession"] == "GCA_000001635.9"])


print("File 'maximal_range.tsv' saved successfully.")
