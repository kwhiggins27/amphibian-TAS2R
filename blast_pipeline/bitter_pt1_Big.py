#!/usr/bin/env python
#SBATCH --job-name=py1B  # Job name
#SBATCH --mail-type=END,FAIL     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kwh1@wi.mit.edu   # Where to send mail
#SBATCH --mem=4gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=py1B_%j.log # Standard output and error log
#SBATCH --error=py1B_%j.err

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

my_path=os.getcwd()
tsv_files=glob.glob(os.path.join(folder+"/blastF/", "*.tsv"))
print(folder)
print(os.path.join(folder+"/blastF/"))
print(tsv_files)

print(len(tsv_files))


desired_colnames = ["name", "gene", "evalue", "s. start", "s. end", "pident","scratch","s. seq", "q. start", "q. end"]
master_list_final=pd.DataFrame(columns=desired_colnames)


mini_list=[]


for f in sorted(tsv_files):
    master_list=pd.read_csv(f, sep = "\t",header = None)
    short_path=os.path.basename(f).rpartition("_local")[0]
    master_list.columns = desired_colnames
    master_list["s. length"] = abs(master_list["s. start"]-master_list["s. end"])
    new_starts=[]
    new_stops=[]

    if borders=="narrow":
        for i, value in master_list.iterrows():
             if value["s. start"] <value["s. end"]:
                counter = counter+1
                new_starts.append(value["s. start"] - (1500- value["s. length"]))
                start_temp = new_starts.append(value["s. start"] - (1500- value["s. length"]))
                new_stops.append(value["s. end"] + (1500- value["s. length"]))
                stop_temp = value["s. end"] + (1500- value["s. length"])
             else:
                counter = counter+1
                new_starts.append(value["s. end"] - (1500- value["s. length"]))
                start_temp = value["s. end"] - (1500- value["s. length"])
                new_stops.append(value["s. start"] + (1500- value["s. length"]))
                stop_temp = value["s. start"] + (1500- value["s. length"])

    elif borders=="new":
        for i, value in master_list.iterrows():
            if value["s. start"] <value["s. end"]:
                if value["q. start"] <= 1:
                    new_starts.append(int(float(value["s. start"])))
                    start_temp = int(float(value["s. start"]))
                    new_stops.append(int(float(value["s. end"])) + (1500- int(float(value["s. length"]))))
                    stop_temp = int(float(value["s. end"])) + (1500- int(float(value["s. length"])))
                else:
                    new_starts.append(int(float(value["s. start"])) - (1500- int(float(value["s. length"]))))
                    start_temp = int(float(value["s. start"])) - (1500- int(float(value["s. length"])))
                    new_stops.append(int(float(value["s. end"])) + (1500- int(float(value["s. length"]))))
            else:
                if value["q. start"] <= 1:
                    new_starts.append(int(float(value["s. end"])) - (1500- int(float(value["s. length"]))))
                    start_temp = int(float(value["s. end"])) - (1500- int(float(value["s. length"])))
                    new_stops.append(int(float(value["s. start"])))
                    stop_temp = int(float(value["s. start"]))
                else:
                    new_starts.append(int(float(value["s. end"])) - (1500- int(float(value["s. length"]))))
                    start_temp = int(float(value["s. end"])) - (1500- int(float(value["s. length"])))
                    new_stops.append(int(float(value["s. start"])) + (1500- int(float(value["s. length"]))))
                    stop_temp = int(float(value["s. start"])) + (1500- int(float(value["s. length"])))
    elif borders=="wide":
        for i, value in master_list.iterrows():
            if value["s. start"] <value["s. end"]:
                new_starts.append(max((value["s. start"])-2000,0))
                start_temp = max((value["s. start"])-2000,0)
                new_stops.append(value["s. end"]+2000)
                stop_temp = value["s. end"]+2000
            else:
                new_stops.append((value["s. start"])+2000)
                start_temp = (value["s. start"])+2000
                new_starts.append(max((value["s. end"])-2000,0))
                stop_temp = max((value["s. end"])-2000,0)
    else:
        print("Problem with borders")

    # master_list["new start"]=int(float(start_temp))
    # master_list["new stop"]=int(float(stop_temp))
    # master_list["new start"]=[int(float(dum)) for dum in new_starts]
    # master_list["new stop"]=[int(float(dum)) for dum in new_stops]
    master_list["new start"]=[int(str(dum).replace(".0", "")) for dum in new_starts]
    master_list["new stop"]=[int(str(dum).replace(".0", "")) for dum in new_stops]
    # master_list["new start"]=pd.Series([int(float(dum)) for dum in new_starts], dtype=int)
    # master_list["new stop"]=pd.Series([int(float(dum)) for dum in new_stops], dtype=int)
    print(master_list.shape)
    # print((master_list["new start"][0]))

    cat=pd.DataFrame({'a':[1,3,4],'b':[4,3,2]})

    #assert False, "%s %s" %(type(master_list),type(master_list["new start"]))
    #assert False, "%s"%(type(cat))

    master_list_final=pd.concat([master_list_final, master_list], axis=0, ignore_index=True)  #ignore_index=True is new
    master_list_final['new start'] = master_list_final['new start'].astype(int)
    master_list_final['new stop'] = master_list_final['new stop'].astype(int)

    header_file = open("/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/pipeline/header_big.txt", "r+")
    header = header_file.read()
    header_file.close

    print("makes it here")

    pull_name= folder + "/pull/" + short_path + "_full_pull.sh"
    genome_chunk_pre = os.path.dirname(genome_location)
    genome_chunk = os.path.join(genome_chunk_pre, f"{accession}_subdivided",f"{short_path}.fasta")
    #genome_chunk = genome_location + "/" + short_path + ".fasta"
    x = open(pull_name, "w+")
    x.write(header)
    x.write ("in="+genome_chunk+"\n")
    x.write ("out="+"expanded_file.fasta"+"\n")
    for i, gene in master_list.iterrows():
        # start_temp = str(gene["new start"]).replace(".0", "")
        # stop_temp = str(gene["new stop"]).replace(".0", "")
        add = "samtools faidx $in " + gene["gene"]+":" + str(gene["new start"]) +str("-") + str(gene["new stop"]) + ">> $out \n"
        #print(add)
        x.write(add)
    x.close()

    os.system("chmod +rwx " + pull_name)

    mini_list.append(short_path)

    #assert False, "%s %s" %(type(master_list),type(master_list["new start"]))

pickle_out=open("master_list.pickle","wb")
pickle.dump(master_list_final, pickle_out)
pickle_out.close()
