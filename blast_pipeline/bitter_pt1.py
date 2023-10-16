#!/usr/bin/env python
#SBATCH --job-name=py1  # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kwh1@wi.mit.edu   # Where to send mail
#SBATCH --mem=4gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=py1_%j.log # Standard output and error log
#SBATCH --error=py1_%j.err


##downloaded at v10 stage

import pandas as pd
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import datetime
from datetime import datetime
import csv
import biolib
import importlib
import pickle
import sys

#Set to temp directory
os.chdir(sys.argv[1])
sys.path.insert(0,(sys.argv[1]))

#Import variables
importlib.import_module('for_py2')
from for_py2 import *

#Set desired wd and read input
os.chdir(folder)


#Define column names for full file (colnames) and for the smaller dataframe that I will create later (desired_colnames)
desired_colnames = ["name", "gene", "evalue", "s. start", "s. end", "pident","scratch","s. seq", "q. start", "q. end"]  #THIS IS THE MAIN ONE

#import data
#master_list = pd.read_csv(BLAST_file,sep = "\t", header = None)
master_list = pd.read_csv("forward_blast.tsv",sep = "\t", header = None)
master_list.columns = desired_colnames

len(master_list)
#master_list[0:5]

# #define a new column for the length of the subject sequence.  This should be much shorter than a gene
master_list["s. length"] = abs(master_list["s. start"]-master_list["s. end"])

# new_start=[]
# new_stop=[]

new_starts=[]
new_stops=[]


if borders=="narrow":
    new_dims = pd.DataFrame(columns = ["new start", "new stop", "new length"],index = range(0,len(master_list.index)))

    counter = -1
    for i, value in master_list.iterrows():
         if value["s. start"] <value["s. end"]:
            counter = counter+1
            new_starts.append(value["s. start"] - (1500- value["s. length"]))
            new_stops.append(value["s. end"] + (1500- value["s. length"]))
         else:
            counter = counter+1
            new_starts.append(value["s. end"] - (1500- value["s. length"]))
            new_stops.append(value["s. start"] + (1500- value["s. length"]))

        # for count, value in master_list.iterrows():
        #      if value["s. start"] <value["s. end"]:
        #         counter = counter+1
        #         new_starts["new start"][counter] =  value["s. start"] - (1500- value["s. length"])
        #         new_starts["new stop"][counter] =  value["s. end"] + (1500- value["s. length"])
        #      else:
        #         counter = counter+1
        #         new_starts["new stop"][counter] =  value["s. start"] + (1500- value["s. length"])
        #         new_starts["new start"][counter] =  value["s. end"] - (1500- value["s. length"])
elif borders=="new":
    for i, value in master_list.iterrows():
        if value["s. start"] <value["s. end"]:
            if value["q. start"] <= 1:
                new_starts.append(value["s. start"])
                new_stops.append(value["s. end"] + (1500- value["s. length"]))
            else:
                new_starts.append(value["s. start"] - (1500- value["s. length"]))
                new_stops.append(value["s. end"] + (1500- value["s. length"]))
        else:
            if value["q. start"] <= 1:
                new_starts.append(value["s. end"] - (1500- value["s. length"]))
                new_stops.append(value["s. start"])
            else:
                new_starts.append(value["s. end"] - (1500- value["s. length"]))
                new_stops.append(value["s. start"] + (1500- value["s. length"]))
elif borders=="wide":
    counter = -1
    for i, value in master_list.iterrows():
        if value["s. start"] <value["s. end"]:
            new_starts.append(max((value["s. start"])-2000,0))
            new_stops.append(value["s. end"]+2000)
        else:
            new_stops.append((value["s. start"])+2000)
            new_starts.append(max((value["s. end"])-2000,0))
else:
    print("Problem with borders")

#new_dims = pd.DataFrame(columns = ["new start", "new stop", "new length"],index = range(0,len(master_list.index)))

# counter = -1
# for count, value in master_list.iterrows():
#      if value["s. start"] <value["s. end"]:
#         counter = counter+1
#         new_starts["new start"][counter] =  value["s. start"] - (1500- value["s. length"])
#         new_starts["new stop"][counter] =  value["s. end"] + (1500- value["s. length"])
#      else:
#         counter = counter+1
#         new_starts["new stop"][counter] =  value["s. start"] + (1500- value["s. length"])
#         new_starts["new start"][counter] =  value["s. end"] - (1500- value["s. length"])


# #add dimensions dataframe to starting dataframe
#master_list = pd.concat([master_list, new_dims],axis = 1)
master_list["new start"]=new_starts
master_list["new stop"]=new_stops
master_list.shape

#master_list[master_list["gene"]=="ctg28876_RHIMB__RM170330.28876"]

#write to fasta
header_file = open("/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/header_pull.txt", "r+")
header = header_file.read()
header_file.close

pickle_out=open("master_list.pickle","wb")
pickle.dump(master_list, pickle_out)
pickle_out.close()

x = open(full_pull, "w+")
x.write(header)
x.write ("in="+genome_location+"\n")
x.write ("out=expanded_file.fasta"+"\n")
#x.write ("out="+expanded_file+"\n")
for i, gene in master_list.iterrows():
    #add = "samtools faidx $in " + gene["gene"] + ">> $out \n"  #.replace(".0", "")
    add = "samtools faidx $in " + gene["gene"]+":" + str(gene["new start"]) +str("-") + str(gene["new stop"]) + ">> $out \n"  #.replace(".0", "")
    #print(add)
    x.write(add)
x.close()

y= open(BLAST_species+"_"+tasted+"_"+"expanded_log_"+run+".txt", "w")
y.write("reference species: "+reference_species+"\n"+"BLAST species: " + BLAST_species+"\n"+tasted+"\n"+run+"\n"+"time written: "+datetime.now().strftime("%Y%m%d%H%M%S"))
y.write(str(master_list.shape)+"\n")
y.close()
#x.write

#os.system("chmod +rwx full_pull12.bat")
