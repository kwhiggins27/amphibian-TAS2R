#!/usr/bin/env python
#SBATCH --job-name=py3  # Job name
#SBATCH --mail-type=NONE    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kwh1@wi.mit.edu   # Where to send mail
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
        query.append(row["name"] + "-" + str(1))
    else:
        dir_num.append("-1")
        query.append(row["name"] + "-" + str(-1))

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

#print(reciprocal_list["query"].head(n=5)+"recip_list")


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

# length_threshold_top = 500
# length_threshold_bottom = 200
# e_threshold=1e-10
# sufficient_length2 = master_list[(master_list["expanded length"]>length_threshold_bottom) & (master_list["expanded length"]< length_threshold_top)]# & (master_list["evalue"]<e_threshold)]

# #sufficient_length_sorted=sufficient_length2.sort_values(by='pident', ascending=False)
# sufficient_length_sorted=sufficient_length2.sort_values(by='evalue', ascending=True)

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

# pickle_in=open("from_cd_forrecip2.pickle","rb")
# from_cd_forrecip=pickle.load(pickle_in)
#
# pickle_in=open("matched.pickle","rb")
# matched=pickle.load(pickle_in)
# if isinstance(from_cd["query"], str):
#     from_cd["query"] = pd.Series([from_cd["query"]])
# mini=pd.Series(from_cd.loc[13562,"query"].replace(".0", "").replace("--1", "").replace("-1", ""))
# if isinstance(matched["query"], str):
#     matched["query"] = pd.Series([matched["query"].replace(".0", "").replace("--1", "").replace("-1", "")])

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



# found = []
# found_sstart=[]
# found_sseq=[]
# for index, row in from_cd_forrecip.iterrows():
#     #found.append(row["query"] in matched["query"].values) #This is the original
#     found.append(
#     row["query"].replace('.0', '').replace('--1', '').replace('-1', '')).isin(
#         matched["query"].replace('.0', '').replace('--1', '').replace('-1', '').values)
# print(matched["query"].head(n=5)+"is matched")
#
# print(from_cd_forrecip["query"].head(n=5)+"is from_cd_forrecip")
#
# from_cd_forrecip["found"] = found
# final_genes = from_cd_forrecip[from_cd_forrecip["found"]==True]
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
# matched_query=[]
# for index, row in matched.iterrows():
#     temp2=row["query"].replace("--1--1","--1").replace("-1-1","-1")
#     matched_query.append(temp2)
# matched=matched.drop(["query"],axis=1)
# matched["query"]=matched_query
#
# reciprocal_query=[]
# for index, row in reciprocal.iterrows():
#     temp2=row["query"].replace("--1--1","--1").replace("-1-1","-1")
#     reciprocal_query.append(temp2)
# reciprocal=reciprocal.drop(["query"],axis=1)
# reciprocal["query"]=reciprocal_query
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
#
# flag=matched[matched["s. seq"].str.startswith("M")]
# flag=matched[matched["q. start"] < 5]
#
# long = []
# for index, row in from_cd.iterrows():
#     long.append(row["query"] in flag["query"].values)
# from_cd["long"]=long
#
# from_cd.to_csv("from_cd_py3.csv")
#
# flagged1=from_cd[from_cd["length"]>1050]
# flagged=flagged1[flagged1["long"]==True]
#
# flagged.to_csv("flagged_py3.csv")
# #Define variable which is subset of from_cd not in flagged
# good1=from_cd[from_cd["length"]<=1050]
# good2=flagged1[flagged1["long"]==False]
# good=pd.concat([good1,good2])
#
# good.to_csv("good.csv")
#
# subset = []
# for index, row in flag.iterrows():
#     subset.append(row["query"] in flagged["query"].values)
# flag["flagged"]=subset
# flag["length"] = abs(flag["s. start"]-flag["s. end"])
#
# matched2=flag[flag["flagged"]==True]
# #matched2=matched2[matched2["evalue"]<1.0E-03]
# matched2=matched2[matched2["length"]>700]
#
# matched2["intended_start"]=0
# matched2["intended_stop"]=0
# matched2["change_start"]=0
#
# for i, row in matched2.iterrows():
#     for j, coordinates in coordinates_bitter.iterrows():
#         if coordinates["gene"]==row["assignment"]:
#             if row["direction"]=="+":
#                 matched2["intended_start"][i]=coordinates["start"]
#                 matched2["intended_stop"][i]=coordinates["stop"]
#                 matched2["change_start"][i]=abs(row["s. start"]-coordinates["start"])
#             else:
#                 matched2["intended_start"][i]=coordinates["stop"]
#                 matched2["intended_stop"][i]=coordinates["start"]
#                 matched2["change_start"][i]=abs(row["s. start"]-coordinates["stop"])
#         else: pass
#     else:
#         pass
#
# flagged["corrected_first"]=0
# flagged["corrected_last"]=0
# flagged["corrected_seqs"]=0
# flagged["corrected_query"]=0
# flagged["corrected_size"]=0
#
# for i, candidate in flagged.iterrows():
#     for j, reciprocal in matched2.iterrows():
#         if candidate["query"]==reciprocal["query"]:
#             if reciprocal["direction"]=="+":
#                 flagged["corrected_first"][i]=int(float(candidate["first"]))+(reciprocal["change_start"])
#                 flagged["corrected_last"][i]=int(float(candidate["last_final"]))
#                 aa_int=int((reciprocal["change_start"])/3)
#                 aa_length=len(candidate["seqs"])
#                 flagged["corrected_seqs"][i]=candidate["seqs"][aa_int:aa_length]
#                 flagged["corrected_query"][i]=BLAST_species+"_"+candidate["chr"]+":"+str(flagged["corrected_first"][i])+"-"+str(flagged["corrected_last"][i])+"-"+str(candidate["direction"])
#                 flagged["corrected_size"][i]=flagged["corrected_last"][i]-flagged["corrected_first"][i]
#             else:
#                 flagged["corrected_first"][i]=int(float(candidate["first"]))
#                 flagged["corrected_last"][i]=int(float(candidate["last_final"]))-(float(reciprocal["change_start"]))
#                 aa_int=int((reciprocal["change_start"])/3)
#                 aa_length=len(candidate["seqs"])
#                 flagged["corrected_seqs"][i]=candidate["seqs"][aa_int:aa_length]
#                 flagged["corrected_query"][i]=BLAST_species+"_"+candidate["chr"]+":"+str(flagged["corrected_first"][i])+"-"+str(flagged["corrected_last"][i])+"-"+str(candidate["direction"])
#                 flagged["corrected_size"][i]=flagged["corrected_last"][i]-flagged["corrected_first"][i]
#         else: pass
#     else:
#         pass
#
# flagged.to_csv("flagged_py3.csv")
#
# flagged2=flagged[["names", "chr", "corrected_first", "last_final", "corrected_seqs","direction", "frame", "strand", "attribute","corrected_query"]]
# flagged2=flagged2.rename({"corrected_first":"first", "corrected_seqs":"seqs", "corrected_query":"query"}, axis="columns")

#corrected_final_list
#corrected_final=pd.concat([good,flagged2])
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

sample_fasta=folder+"/for_TMHMM.fasta"
sample_pred=folder+"/TMbed.pred"

os.chdir("/lab/solexa_weng/playground/Kate_Higgins/tmbed/tmbed/")
# os.system("python -m /lab/solexa_weng/playground/Kate_Higgins/tmbed/tmbed/tmbed predict -f %s -p %s"%(sample_fasta, sample_pred))
#                     /lab/solexa_weng/playground/Kate_Higgins/tmbed/tmbed/tmbed
os.system("python -m tmbed predict -f %s -p %s"%(sample_fasta, sample_pred))

os.chdir(sys.argv[1])

file_path=folder+"/TMbed.pred"
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


# os.system("grep '^# ' " + tmhmm_out + "> " +  tmhmm_1)
# os.system("sed '/Length:/d' ./" + tmhmm_1 + " > " +  tmhmm_final)
#
# #sed '/Length:/d' ./TMRs.txt > TMRs3.txt
#
# TMRs = []
# with open(tmhmm_final ) as f:
#     lines = f.readlines()
#
# for i in lines:
#     TMRs.append(i[-2])
#
# TMRs[0:5]
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

# found2 = []
# for index, row in from_cd_forrecip.iterrows():
#     found2.append(row["query"] in final_genes_7TM["query"].values)

# last_final=[]
# direction_final=[]
#
# for index, row in from_cd_forrecip.iterrows():
#     stop=row["last"].find("-", 0)
#     temp=row["last"][0:stop]
#     end_string=len(row["last"])
#     last_final.append(temp)
#     if row["last"][end_string-3:end_string]=="--1":
#         direction_final.append("-")
#     else:
#         direction_final.append("+")
#
# from_cd_forrecip["found2"]=found2
# from_cd_forrecip["query_short"]=0
# from_cd_forrecip["assignment"]=0
# from_cd_forrecip["direction"]=direction_final
# from_cd_forrecip["last_final"]=last_final
# from_cd_forrecip["source"] = "Genome"
# from_cd_forrecip["feature"] = "exon"
# from_cd_forrecip["score"] = "."
#
# final_for_gtf = from_cd_forrecip[from_cd_forrecip["found2"]==True]



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

# # Open the file in append mode
# with open(logfile, 'a', newline='') as csvfile:
#     writer = csv.writer(csvfile)
# 
#     # Write the variables to a new row
#     writer.writerow([accession, latin, common, taxa, number_tas2rs])

    # echo "accession=\"$2\"">> $for_py
    # echo "latin=\"$3\"">> $for_py
    # echo "common=\"$4\"">> $for_py
    # echo "taxa=\"$5\"">> $for_py
    # echo "folder=\"/lab/wengpj01/vertebrate_pipeline/subdirs/${shortname}/$now\"" >> $for_py
    # echo "logfile=\"/lab/wengpj01/vertebrate_pipeline/20230526_summaryruns\""
