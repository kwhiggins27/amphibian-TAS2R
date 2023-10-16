#!/usr/bin/env python
#SBATCH --job-name=py_exp  # Job name
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kwh1@wi.mit.edu   # Where to send mail
#SBATCH --mem=4gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=logs/py_exp_%j.log # Standard output and error log
#SBATCH --error=logs/py_exp_%j.err

import pandas as pd
import os
from os.path import exists
# import glob
# from Bio import SeqIO
# from Bio.Seq import Seq
# import numpy as np
from functools import reduce
import sys
import argparse


species = sys.argv[1]
directory = sys.argv[2]
BLAST_file = sys.argv[3]
full_file=sys.argv[4]

blast = pd.DataFrame()
blast = pd.read_csv(BLAST_file)

add_me = []
for index,row in blast.iterrows():
    add_me.append(row["query"])
blast["Unique"] = add_me

addresses={}



for x in range(1,22):
    if species=="xenopus":
        base=21
        addresses[x]=(f"{429+x}_S{base+x}_L002_fc.tsv")
    if species=="terribilis":
        base=0
        addresses[x]=(f"{408+x}_S{base+x}_L001_fc.tsv")
    if species=="cane":
        base=42
        addresses[x]=(f"{450+x}_S{base+x}_L003_fc.tsv")
    if species=="bullfrog":
        base=63
        addresses[x]=(f"{471+x}_S{base+x}_L004_fc.tsv")
    if species=="axolotl":
        base=9
        addresses[x]=(f"{510+x}_S{base+x}_L004_fc.tsv")

print (addresses)

os.chdir(directory)

xen1 = {}
xen2 = {}

for x in list(range(1,len(addresses)+1)):
    xen1[x]=pd.read_csv(addresses[x], sep="\t",skiprows=1)#[["/lab/wengpj01/xenopus/star_expression_xenopus/430_S22_L002Aligned.sortedByCoord.out.bam"]]

    count = 0
    dummy_col = [False for dum in range(len(xen1[x]))]
    i=0
    for idx,row_fc in xen1[x].iterrows():
        j=0
        for jdx,row_blast in blast.iterrows():
            var1 = int(row_blast["first"])-500
            var2 = int(row_fc["Start"])
            var3 = int(row_fc["End"])
            var4 = int(row_blast["last_final"])+500
            var5 = row_blast["chr"]
            var6 = row_fc["Chr"]
            if var5==var6:
                if var1 <= var2 and var1 <= var3 and var4 >= var2 and var4 >= var3:
                    dummy_col[i]=True
                    count = count + 1
                    j+=1
                    continue

            else:
                pass
            j+=1
        i+=1

        fc2=xen1[x]
        xen1[x]["in_blast"] = dummy_col


xen3 = {}
for x in list(range(1,len(addresses)+1)):
    xen2[x] = xen1[x][xen1[x]["in_blast"]==True]
    unique = []
    for y, z in xen2[x].iterrows():
        unique.append(z["Chr"] + ":" + str(z["Start"]) + "-" + str(z["End"]))
    xen2[x]["Unique"]=unique
    long_name = xen2[x].keys()[6]
    xen3[x] = xen2[x].rename({xen2[x].keys()[6]: "Counts"}, axis="columns")
    xen3[x]=xen3[x][["Unique", "Counts"]]


xen3[1].columns=["id", "Counts_1.1"]
xen3[2].columns=["id", "Counts_1.2"]
xen3[3].columns=["id", "Counts_1.3"]
xen3[4].columns=["id", "Counts_1.4"]
xen3[5].columns=["id", "Counts_1.5"]
xen3[6].columns=["id", "Counts_1.6"]
xen3[7].columns=["id", "Counts_1.7"]
xen3[8].columns=["id", "Counts_2.1"]
xen3[9].columns=["id", "Counts_2.2"]
xen3[10].columns=["id", "Counts_2.3"]
xen3[11].columns=["id", "Counts_2.4"]
xen3[12].columns=["id", "Counts_2.5"]
xen3[13].columns=["id", "Counts_2.6"]
xen3[14].columns=["id", "Counts_2.7"]
xen3[15].columns=["id", "Counts_3.1"]
xen3[16].columns=["id", "Counts_3.2"]
xen3[17].columns=["id", "Counts_3.3"]
xen3[18].columns=["id", "Counts_3.4"]
xen3[19].columns=["id", "Counts_3.5"]
xen3[20].columns=["id", "Counts_3.6"]
xen3[21].columns=["id", "Counts_3.7"]



short_names=[]
for x in range(1,22):
    #print(f"{272+x}_S{base+x}")
    short_names.append((f"xen1[{x}]"))


xen4=reduce(lambda  left,right: pd.merge(left,right,on=['id'],how='outer'), [xen3[k] for k in sorted(xen3.keys())])


xen4.to_csv(full_file)
