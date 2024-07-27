#!/usr/bin/env python
#SBATCH --job-name=py3  # Job name
#SBATCH --mail-type=NONE    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremailaddress@yourinstitute  # Where to send mail
#SBATCH --mem=100gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=py3_%j.log # Standard output and error log
#SBATCH --error=py3_%j.err

print("checkpoint -1")

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

print("checkpoint 0")

#Set to temp directory
os.chdir(sys.argv[1])
sys.path.insert(0, (sys.argv[1]))

#Import variables
importlib.import_module('for_py2')
from for_py2 import *

#Set desired wd and read input
os.chdir(folder)

#import pickle
pickle_in=open("from_cd_forrecip.pickle","rb")
from_cd_forrecip=pickle.load(pickle_in)

from_cd_forrecip.to_csv("from_cd_forrecip_check.csv")
# ## Reverse blast

coordinates_bitter=pd.read_csv(coordinates,header=0, names=["gene", "chromosome", "start", "stop", "scratch"], )
with open(coordinates, 'r') as file:
    data = file.read()

# #Define column names for full file (colnames) and for the smaller dataframe that I will create later (desired_colnames)
#desired_colnames_reciprocal = ["name", "gene1", "evalue2", "s. start2", "s. end2", "pident","scratch","s. seq2"]
#desired_colnames = ["name", "gene", "evalue", "s. start", "s. end", "pident","s. seq"]
desired_colnames = ["name", "gene", "evalue", "s. start", "s. end", "pident","scratch","s. seq", "q. start", "q. end"]
#desired_colnames = ["name", "gene", "evalue", "q. start", "q. end", "s. start", "s. end", "s. seq"]

#import data
reciprocal_list = pd.read_csv("reciprocal_blast.tsv", sep = "\t",header = None)
reciprocal_list.columns = desired_colnames

reciprocal_list["query_short"]=0
reciprocal_list["assignment"]=0
reciprocal_list["direction"]=0

dir_num=[]
query=[]
for index, row in reciprocal_list.iterrows():
    if row["s. start"] < row["s. end"]:
        dir_num.append("1")
        query.append(str(row["name"]) + str("-") + str(1))
    else:
        dir_num.append("-1")
        query.append(str(row["name"]) + str("-") + str(-1))

reciprocal_list["dir_num"]=dir_num
reciprocal_list["query"]=query

print(query[0:20])


for i, gene in enumerate(reciprocal_list["name"]):
    if reciprocal_list["s. start"][i] < reciprocal_list["s. end"][i]:
        reciprocal_list["direction"][i]="+"
        # reciprocal_list["dir_num"][i]="1"
    else:
        reciprocal_list["direction"][i]="-"
        # reciprocal_list["dir_num"][i]="-1"

query = []

print(reciprocal_list["query"][:20])


#Which reciprocal_query fall within coordinates?
margin = 500
inner_margin = 5
for i, gene in reciprocal_list.iterrows():
    for j, coordinates in coordinates_bitter.iterrows():
        if gene["gene"]== coordinates["chromosome"]:
            if ((((int(float(gene["s. start"])) + inner_margin) > (int(float(coordinates["start"]))-margin)) &  ((int(float(gene["s. start"])) -margin) < (int(float(coordinates["stop"]))+inner_margin))) or ((inner_margin+(int(float(gene["s. end"])))> (int(float(coordinates["start"]))-margin)) &  ((int(float(gene["s. end"]))-margin)< (int(float(coordinates["stop"]))+ inner_margin)))):
                reciprocal_list["assignment"][i]=coordinates["gene"]
            else: pass
        else:
            pass

matched = reciprocal_list[reciprocal_list["assignment"]!=0]
# matched['query'] = matched['query']

pickle_out0=open("reciprocal_list.pickle","wb")
pickle.dump(reciprocal_list, pickle_out0)
pickle_out0.close()

pickle_in=open("reciprocal_list.pickle","rb")
reciprocal_list=pickle.load(pickle_in)

pickle_out0=open("matched1.pickle","wb")
pickle.dump(matched, pickle_out0)
pickle_out0.close()

# #write to fasta
x = open(BLAST_species+"_reciprocal_expanded.fasta", "w")
for i, gene in matched.iterrows():
    add = ">"+str(gene["query"]) + "\n" +str(gene["s. seq"]).replace("-","") + "\n"
    #print(add)
    x.write(add)
x.close()

#remove exact duplicates
threshold = str(1)
input_file = BLAST_species+"_reciprocal_expanded.fasta"
cdhit_output = BLAST_species+"_reciprocal_cdhit_"+threshold+".fasta"


os.system("cd-hit -i " + input_file + " -o " + cdhit_output + " -c "+ threshold)

#read in output file as "expanded."  Note that the order of entries has been preserved
cd_out_seqs_r = []
for seq in SeqIO.parse(cdhit_output, "fasta"):
    cd_out_seqs_r.append(str(seq.seq))
cd_out_names_r=[]
cd_out_file_r=open(cdhit_output, "r")
#print(cd_out_file.readlines())
for i in range(0,len(cd_out_seqs_r)):
    cd_out_file=open(cdhit_output, "r")
    cd_out_names_r.append(cd_out_file.readlines()[i*2])

print(cd_out_names_r[0:5])

cd_out_names_short_r=[]
for i, name in enumerate(cd_out_names_r):
    start=name.find("_",0)
    stop=name.find("_", start+1)
    cd_out_names_short_r.append(name[1:stop])

print(cd_out_names_short_r[0:5])

matched = reciprocal_list[reciprocal_list["assignment"]!=0]
print(matched.shape)

print(from_cd_forrecip["query"][:20])



pickle_out1=open("from_cd_forrecip2.pickle","wb")
pickle.dump(from_cd_forrecip, pickle_out1)
pickle_out1.close()

pickle_out2=open("matched2.pickle","wb")
pickle.dump(matched, pickle_out2)
pickle_out2.close()

from_cd_query_split = from_cd_forrecip["query"].str.split(":", n=1, expand=True)
from_cd_query_split[1] = from_cd_query_split[1].str.replace(".0", "").str.replace("--1", "").str.replace("-1", "")
from_cd_query = from_cd_query_split[0] + ":" + from_cd_query_split[1]
from_cd_query = from_cd_query.reset_index(drop=True)
matched_query_split = matched["query"].str.split(":", n=1, expand=True)
matched_query_split[1] = matched_query_split[1].str.replace(".0", "").str.replace("--1", "").str.replace("-1", "")
matched_query = matched_query_split[0] + ":" + matched_query_split[1]
matched_query= matched_query.reset_index(drop=True)

print(from_cd_query[0])
print(matched_query[0])

subset = from_cd_query.isin(matched_query)
subset=subset.reset_index(drop=True)
from_cd_forrecip["found"] = subset.values

final_genes = from_cd_forrecip[from_cd_forrecip["found"] == True]

print(final_genes.shape)
final_genes.shape
lost_genes = from_cd_forrecip[from_cd_forrecip["found"]==False]
print(lost_genes.shape)

x = open(BLAST_species+"_"+tasted+"_"+"confirmed_protein_"+run+".fasta", "w")
for i, gene in final_genes.iterrows():
    add =  ">"+BLAST_species+"_"+str(gene["chr"])+":"+str(gene["first"])+"-"+str(gene["last"])+"\n" +str(gene["seqs"]).replace("-","") + "\n"
    #print(add)
    x.write(add)
x.close()

y = open(BLAST_species+"_"+tasted+"_"+"lost_protein_"+run+".fasta", "w")
for i, gene in lost_genes.iterrows():
    add =  ">"+BLAST_species+"_"+str(gene["chr"])+":"+str(gene["first"])+"-"+str(gene["last"])+"\n" +str(gene["seqs"]).replace("-","") + "\n"
    #print(add)
    y.write(add)
y.close()

final_genes.to_csv("passed_blastR.csv")

# ## Deal with edge cases
from_cd=final_genes
reciprocal=reciprocal_list
#
last_final=[]
length=[]
query_fixed=[]

for index, row in from_cd.iterrows():
    stop=row["last"].find("-", 0)
    temp=row["last"][0:stop]
    temp2=row["query"].replace("--1--1","--1").replace("-1-1","-1").replace(".0","")
    last_final.append(temp)
    length.append(int(float(temp))- int(float(row["first"])))
    query_fixed.append(temp2)
from_cd["last_final"]=last_final
from_cd["length"]=length
from_cd=from_cd.drop(["query"],axis=1)
from_cd["query"]=query_fixed

corrected_final=final_genes
x = open("for_TMHMM.fasta", "w")
for i, gene in corrected_final.iterrows():
    add =  ">"+str(gene["query"])+"\n" +str(gene["seqs"]).replace("-","") + "\n"
    #print(add)
    x.write(add)
x.close()


pickle_out=open("pre_TMHMM.pickle","wb")
pickle.dump(corrected_final, pickle_out)
pickle_out.close()


#Run TMbed

sample_fasta="../../../subdirs" + accession + "/for_TMHMM.fasta"
sample_pred="../../../"+accession+"/TMbed.pred"

os.chdir("tmbed/tmbed/")
os.system("python -m tmbed predict -f %s -p %s"%(sample_fasta, sample_pred))

os.chdir(sys.argv[1])

file_path="TMbed.pred"
with open(file_path, 'r') as file:
    file_lines = file.read().splitlines()

header = []
seq = []
topo = []
for i in range(0, len(file_lines), 3):
    header.append(file_lines[i][1:])
    seq.append(file_lines[i+1])
    topo.append(file_lines[i+2].lower())

data = {'header': header, 'seq': seq, 'topo': topo}
df = pd.DataFrame(data)

# Count the occurrences of "..hh" in each row and save as a new column 'TMs'
df['TMs'] = df['topo'].str.count(r'\.\.hh')

TMRs = df['TMs'].astype(str).tolist()
#
corrected_final["TMRs"] = TMRs

corrected_final.to_csv("post_TMbed.csv")

final_genes_7TM = corrected_final[corrected_final["TMRs"]==str(7)]
lost2 = corrected_final[corrected_final["TMRs"]!=str(7)]
lost2.shape

final_genes_7TM.to_csv("final_genes_7TM.csv")

x = open(BLAST_species+"_"+tasted+"_"+"protein_confirmed2x"+run+".fasta", "w")
for i, gene in final_genes_7TM.iterrows():
    add =  ">"+BLAST_species+"_"+str(gene["chr"])+":"+str(gene["first"])+"-"+str(gene["last_final"])+"\n" +str(gene["seqs"]).replace("-","") + "\n"
    #print(add)
    x.write(add)
x.close()

y = open(BLAST_species+"_"+tasted+"_"+"not7TM_"+run+".fasta", "w")
for i, gene in lost2.iterrows():
    add =  ">"+BLAST_species+"_"+str(gene["chr"])+":"+str(gene["first"])+"-"+str(gene["last_final"])+"\n" +str(gene["seqs"]).replace("-","") + "\n"
    #print(add)
    y.write(add)
y.close()

pickle_out=open("final_genes_7M.pickle","wb")
pickle.dump(final_genes_7TM, pickle_out)
pickle_out.close()

final_genes_7TM["source"]="Genome"
final_genes_7TM["feature"]="exon"
final_genes_7TM["score"]="."

#
final_for_gtf_rearranged = final_genes_7TM[["chr", "source", "feature", "first", "last_final", "score", "strand", "frame", "attribute"]]
final_for_gtf_rearranged.to_csv (BLAST_species + "_from_pipeline.gtf", sep="\t", header=False, index=False,quoting=csv.QUOTE_NONE)

length = []
for index, rows in final_for_gtf_rearranged.iterrows():
    length.append(int(float(rows["last_final"])) - int(float(rows["first"])))

number_tas2rs=len(final_for_gtf_rearranged.index)

# Open the file in append mode
with open(logfile, 'a', newline='') as csvfile:
    writer = csv.writer(csvfile)

    # Write the variables to a new row
    writer.writerow([accession, latin, common, taxa, number_tas2rs])
