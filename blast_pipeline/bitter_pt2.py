#!/usr/bin/env python
#SBATCH --job-name=py2  # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kwh1@wi.mit.edu   # Where to send mail
#SBATCH --mem=200gb          # Job memory request, down from 200 and closer to the 64gb I think you're using per instance
#SBATCH --nodes=1           # ensure cores are on one node
#SBATCH --ntasks=1          # run a single task
#SBATCH --cpus-per-task=1       # number of cores/threads requested, up from 4 you're asking for now - see too that I changed --runThreadN below to match
#SBATCH --partition=20
#SBATCH --output=py2_%j.log # Standard output and error log
#SBATCH --error=py2_%j.err

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

#import pickle
pickle_in=open("master_list.pickle","rb")
master_list=pickle.load(pickle_in)

print(str(master_list.shape) + " is shape of master_list from pickle")


# #read in output file as "expanded."  Note that the order of entries has been preserved
expanded_file="expanded_file.fasta"
expanded = []
for seq in SeqIO.parse(expanded_file, "fasta"):
#for seq in SeqIO.parse(expanded_file, "fasta"):
    expanded.append(seq.seq)

print(str(len(expanded))+ " is length of expanded")

#expanded

#Massive loop that steps through each sequence in all possible reading frames and identifies the longest possible reading frame.  Output is all_long_seqs and all_long_lengths.

#for each sequence in expanded
all_long_seqs = []
all_long_lengths = []
all_starts=[]
all_stops=[]
all_lengths_check=[]
direction=[]
k_count=[]
longest_read_na_seq=[]
nt_length=[]
all_k_len=[]
seq01=[]
trans0=[]
for i, seq in enumerate(expanded):
    #Define 6 reading frames
    seq0=Seq(expanded[i]) #this is new
    seq1 = Seq("TGATGA" + seq0)
    seq2 = Seq("TGATGAN" + str(seq0))
    seq3 = Seq("TGATGANN" + str(seq0))
    seq4 = Seq("TGATGA" + seq0.reverse_complement())
    seq5 = Seq("TGATGAN" + str(seq0.reverse_complement()))
    seq6 = Seq("TGATGANN" + str(seq0.reverse_complement()))
    seqs = [seq1, seq2, seq3, seq4, seq5, seq6]

    #translate them
    trans = []
    for j,na in enumerate(seqs):
        trans.append(str(na.translate())+"*")

    seq01.append(str(seq2)[0:30])
    trans0.append(str(trans[1])[0:30])

    longest_read_length = 0
    longest_read_seq = ""
    longest_read_count = 0
    longest_start=0
    longest_stop=0
    direction_temp=0
    k_temp=0

    #main loop
    longest_reading_frame = ""
    for k, aa in enumerate(trans):
        current_loc = 0
        #loop across all stop codons
        counter = 0

        for translated in range(aa.count("*")):
            ##counter which stop codon
            counter+= 1
            current_loc = aa.find("*",current_loc+1, len(aa))
            next_M = aa.find("M", current_loc, len(aa))
            next_stop = aa.find("*", current_loc+1, len(aa)-1)
            current_length = next_stop - next_M
        #     print(current_loc)
        #     print(next_M)
        #     print(str(next_stop) +"\n")
            if next_M>next_stop:
                pass
            elif next_M == -1:
                pass
            else:
                if current_length>longest_read_length:
                    longest_read_length = current_length
                    longest_read_count = counter
                    longest_read_seq = aa[next_M:next_stop+1]
                    longest_reading_frame = k
                    #longest_start=next_M
                    #longest_stop=next_stop
                    #na_seq=""
                    k_temp=k

                    if k<3:
                        direction_temp=1
                    else:
                        direction_temp=-1
                    if k==0:
                        longest_start=next_M*3 +2
                        longest_stop=next_stop*3 +2
                        na_seq = str(seq1[longest_start:longest_stop])
                    elif k==1:
                        longest_start=next_M*3 +1  #was +1
                        longest_stop=next_stop*3 +1
                        na_seq = str(seq2[longest_start:longest_stop])
                    elif k==2:
                        longest_start=next_M*3  #was 0
                        longest_stop=next_stop*3
                        na_seq = str(seq3[longest_start:longest_stop])
                    elif k==3:
                        longest_start=next_M*3 +2
                        longest_stop=next_stop*3 +2
                        na_seq = str(seq4[longest_start:longest_stop])
                    elif k==4:
                        longest_start=next_M*3 +1
                        longest_stop=next_stop*3 +1
                        na_seq = str(seq5[longest_start:longest_stop])
                    elif k==5:
                        longest_start=next_M*3
                        longest_stop=next_stop*3
                        na_seq = str(seq6[longest_start:longest_stop])
                    else:
                        longest_start=-100
                        longest_stop=-100
    all_long_seqs.append(longest_read_seq)
    all_long_lengths.append(longest_read_length)
    all_starts.append(longest_start-6)
    all_stops.append(longest_stop-6)
    all_lengths_check.append(longest_stop-longest_start)
    direction.append(direction_temp)
    k_count.append(k_temp)
    longest_read_na_seq.append(na_seq)
    nt_length.append(len(seq1))

#add output back to starting dataframe
master_list["expanded seq"] = all_long_seqs
master_list["expanded length"] = all_long_lengths
master_list["direction"]= direction
master_list["offset start"]=all_starts
master_list["offset stop"]=all_stops
# master_list["true start"]= 0
# master_list["true stop"]=0
# master_list["true length"]=0
master_list["k"]=k_count
master_list["na_seq"]=longest_read_na_seq
master_list["nt_length"]=nt_length
master_list["trans0"]=trans0
master_list["seq01"] = seq01


first=[]
last=[]
truncation=[]
t_left=[]
t_right=[]
for index, row in master_list.iterrows():
    if row["s. start"] <row["s. end"]:
        req_length=int(float(row["new start"]))-int(float(row["new stop"]))
        trunc=row["nt_length"]+(int(float(row["new start"])))-(int(float(row["new stop"])))
        trunc_left= abs(min(0, row["new start"]-1))
        t_left.append(trunc_left)
        trunc_right=max(req_length - row["nt_length"] - trunc_left,0)
        t_right.append(trunc_right)
        a = row["new start"] + trunc_left + row["offset start"]-2
        b = a + (3*row["expanded length"]+3) -1
#         a = max((int(float(row["new start"])))+int(float(row["offset start"]))-2,1)
#         b = max((int(float(row["new start"])))+int(float(row["offset stop"])),1)
        first.append(a)
        last.append(b)

    else:
        req_length=int(float(row["new stop"]))-int(float(row["new start"]))
        trunc=(int(float(row["new stop"])))-(int(float(row["new start"])))-row["nt_length"]
        trunc_left= abs(min(0, row["new start"]-1))  #had new start here previously, start of Dec
        t_left.append(trunc_left)
        trunc_right=max(req_length - row["nt_length"] +7 - trunc_left,0)
        t_right.append(trunc_right)
        a = row["new stop"] - trunc_right - row["offset start"]+2
        b = a - (3*row["expanded length"]+3) +1

#             a = max((int(float(row["new stop"])))-int(float(row["offset start"]))-trunc+1,1)
#             b = max((int(float(row["new stop"])))-int(float(row["offset stop"]))-trunc-1,1)
        first.append(b)
        last.append(a)

master_list["first"] = first
master_list["last"] = last
master_list["true length"] = master_list["last"] -master_list["first"]
master_list["trunc left"]=t_left
master_list["trunc right"]=t_right

master_list.to_csv("master_list.csv")

# #write to fasta NUCLEOTIDE
length_threshold_bottom = int(min_aa)*3 #150*3
length_threshold_top = int(max_aa)*3
e_threshold=1.0e-10#1.0e-10
sufficient_length2 = master_list[(master_list["true length"]>length_threshold_bottom) & (master_list["true length"]< length_threshold_top) & (master_list["evalue"]<e_threshold)]

print(length_threshold_top)


#remove if has more than 10 consecutive Ns
#sufficient_length_complete=[x for x in master_list["na_seq"] if "NNNNNNNNNN" not in x]
sufficient_length_complete=sufficient_length2[sufficient_length2["na_seq"].str.contains("NNNNNNNNNNNNNNNNNNNN")==False]
sufficient_length_complete2=sufficient_length_complete[sufficient_length_complete["na_seq"].str.contains("nnnnnnnnnnnnnnnnnnnn")==False]

#sufficient_length_sorted=sufficient_length2.sort_values(by='pident', ascending=False)
sufficient_length_sorted=sufficient_length_complete2.sort_values(by='evalue', ascending=True)

#sufficient_length_sorted.shape

x = open(BLAST_species+"_"+tasted+"_"+"expanded_nucleotide_"+run+".fasta", "w")
for i, gene in sufficient_length_sorted.iterrows():
    add = ">"+BLAST_species+"_"+str(gene["gene"])+":"+str(gene["first"])+"-"+str(gene["last"])+"-"+str(gene["direction"])+"\n" +str(gene["na_seq"]).replace("-","") + "\n"
    #print(add)
    x.write(add)
x.close()

y= open("log_"+BLAST_species+"_"+tasted+"_"+"expanded_nucleotide_"+run+".txt", "w")
y.write("reference species: "+reference_species +"\n"+"BLAST species: " + BLAST_species +"\n"+tasted+"\n"+run+"\n"+"time written: "+datetime.now().strftime("%Y%m%d%H%M%S"))
y.write(str(master_list.shape)+"\n")
y.write(expanded_file +"\n")
y.write(str(sufficient_length_sorted.shape))
y.close()

#master_list["s. seq"][2]
sufficient_length_sorted.shape

#remove exact duplicates
threshold = str(1.0)
input_file = BLAST_species+"_"+tasted+"_"+"expanded_nucleotide_"+run+".fasta"
cdhit_output = BLAST_species+"_"+tasted+"_"+"expanded_nucleotide_"+run+threshold+".fasta"


os.system("cd-hit-est -i " + input_file + " -o " + cdhit_output + " -c "+ threshold)

##### #write to fasta PEPTIDE

length_threshold_bottom = int(min_aa)*3 #150*3
length_threshold_top = int(max_aa)*3
e_threshold=1.0e-10
sufficient_length2 = master_list[(master_list["true length"]>length_threshold_bottom) & (master_list["true length"]< length_threshold_top) & (master_list["evalue"]<e_threshold)]
sufficient_length2=sufficient_length2[sufficient_length2["expanded seq"].str.contains("XXXXXXXXXXXXXXXXXXXX")==False]
sufficient_length2=sufficient_length2[sufficient_length2["expanded seq"].str.contains("xxxxxxxxxxxxxxxxxxxx")==False]

#sufficient_length_sorted=sufficient_length2.sort_values(by='pident', ascending=False)
sufficient_length_sorted=sufficient_length2.sort_values(by='evalue', ascending=True)

sufficient_length_sorted.shape
time2=datetime.now().strftime("%Y%m%d%H%M%S")
x = open(BLAST_species+"_"+tasted+"_"+"expanded_protein_"+run+".fasta", "w")
for i, gene in sufficient_length_sorted.iterrows():
    add =  ">"+BLAST_species+"_"+str(gene["gene"])+":"+str(gene["first"])+"-"+str(gene["last"])+"-"+str(gene["direction"])+"\n" +str(gene["expanded seq"]).replace("-","") + "\n"
    #print(add)
    x.write(add)
x.close()

len(sufficient_length_sorted)

# search = sufficient_length_sorted[(sufficient_length_sorted["gene"]=="NC_030685.2") & (sufficient_length_sorted["first"]<12031400) & (sufficient_length_sorted["first"]>12031300)]
# search = search.reset_index()# &

#remove exact duplicates
threshold = str(1.0)
input_file = BLAST_species+"_"+tasted+"_"+"expanded_protein_"+run+".fasta"
cdhit_output = "for_blastR.fasta"


try:
    command = "cd-hit -i " + input_file + " -o " + cdhit_output + " -c " + str(threshold)
    os.system(command)
except Exception as e:
    # print("Command execution failed. Writing DataFrame to CSV instead.")
    with open(logfile, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the variables to a new row
        writer.writerow([accession, latin, common, taxa, "0"])
        with open("/lab/wengpj01/vertebrate_pipeline/list_complete.txt", "a") as file:
            file.write(accession + "\n")
        sys.exit(1)

try:
    assert os.path.exists(cdhit_output), "The file does not exist."
    print("The file exists.")
except AssertionError as e:
        with open(logfile, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Write the variables to a new row
            writer.writerow([accession, latin, common, taxa, "0"])
            with open("/lab/wengpj01/vertebrate_pipeline/list_complete.txt", "a") as file:
                file.write(accession + "\n")
            sys.exit(1)
except Exception as e:
    print("An error occurred while checking the file existence:", str(e))
    sys.exit(1)


# os.system("cd-hit -i " + input_file + " -o " + cdhit_output + " -c "+ threshold)

y= open("log_"+BLAST_species+"_"+tasted+"_"+"expanded_protein_"+run+".txt", "w")
y.write("reference species: "+reference_species +"\n"+"BLAST species: " + BLAST_species +"\n"+tasted+"\n"+run+"\n"+"time written: "+datetime.now().strftime("%Y%m%d%H%M%S"))
y.write(str(master_list.shape)+"\n")
y.write(expanded_file +"\n")
y.write(str(sufficient_length_sorted.shape))
y.close()

print("Success")
##Temp redefine cdhit_output to match reciprocal blast
#cdhit_output = "/lab/solexa_weng/playground/Kate_Higgins/other_vertebrates/xenopus/xenopus_bitter_expanded_protein_v1_20230401105724_1.0.fasta"


#read in output file as "expanded."  Note that the order of entries has been preserved
cd_out_seqs = []
for seq in SeqIO.parse(cdhit_output, "fasta"):
    cd_out_seqs.append(str(seq.seq))
cd_out_names=[]
cd_out_file=open(cdhit_output, "r")
#print(cd_out_file.readlines())
for i in range(0,len(cd_out_seqs)):
    cd_out_file=open(cdhit_output, "r")
    cd_out_names.append(cd_out_file.readlines()[i*2])

cd_out_names_short=[]
pull_chr = []
pull_start=[]
pull_stop=[]
if chrom_name_have_underscore=="yes":
    for i, name in enumerate(cd_out_names):
        #start=name.find("_",0)
        start1=name.find("_", 1) #after cane
        start2=name.find("_", start1+1) #after chsome
        start3=name.find(":", start2)
        start4=name.find("-", start3)
        start5=name.find("\n", start4)
        #cd_out_names_short.append(name[start+1:stop])
        cd_out_names_short.append(name[1:start1])
        if BLAST_species=="terribililis":
            pull_chr.append(name[start1+1:start2])
        else:
            pull_chr.append(name[start1+1:start3])
        pull_start.append(name[start3+1:start4])
        pull_stop.append(name[start4+1:start5])
elif chrom_name_have_underscore=="no":
    for i, name in enumerate(cd_out_names):
        #start=name.find("_",0)
        start1=name.find("_", 1) #after cane
        start3=name.find(":", start1)
        start4=name.find("-", start3)
        start5=name.find("\n", start4)
        #cd_out_names_short.append(name[start+1:stop])
        cd_out_names_short.append(name[1:start1])
        if BLAST_species=="terribililis":
            pull_chr.append(name[start1+1:start2])
        else:
            pull_chr.append(name[start1+1:start3])
        pull_start.append(name[start3+1:start4])
        pull_stop.append(name[start4+1:start5])
else:
    print("Problem with chromosome nomenclature")

from_cd=pd.DataFrame()
from_cd["names"]=cd_out_names_short
from_cd["chr"] = pull_chr
from_cd["first"] = pull_start
from_cd["last"] = pull_stop
from_cd["seqs"]=cd_out_seqs
from_cd["unique"]=from_cd["names"]+"_"+from_cd["seqs"]
cd_query =[]
for i in cd_out_names:
    cd_query.append(i[1:-1])
from_cd["query"] = cd_query

from_cd.to_csv("from_cd.csv")

# from_cd["query"] = cd_out_names[0][1:-2]

# master_list["unique"]=str(0)
# for i, name in enumerate(master_list["name"]):
#     master_list["unique"][i]=str(str(name).replace("Frog_",""))+"_"+str(master_list["expanded seq"][i])

master_list["seqs"] = master_list["expanded seq"]

pickle_out=open("from_cd_troubleshooting.pickle","wb")
pickle.dump(from_cd, pickle_out)
pickle_out.close()

# from_cd["gene"]=str(0)
# from_cd["new start"]=str(0)
# from_cd["new stop"]=str(0)

from_cd_merged=from_cd.merge(master_list, how="left", on="seqs")
print(len(from_cd_merged))
len(from_cd)

#Note: there are two receptors with the same origin and same sequence right next to eachother
#print(from_cd_merged[21:23])

from_cd_merged = from_cd_merged.sort_values("expanded length", ascending=False).drop_duplicates("seqs").sort_index()

# from_cd_forrecip = from_cd_merged [["names", "chr", "first_x", "last_x", "seqs", "direction", "query"]]
from_cd_merged=from_cd_merged[["names", "chr", "first_x", "last_x", "direction", "seqs"]]
from_cd_merged=from_cd_merged.sort_values(["chr", "first_x"], ascending=(True, True))

from_cd_merged["source"] = "Genome"
from_cd_merged["feature"] = "exon"
from_cd_merged["score"] = "."
from_cd_merged["strand"] = "0"
strand =[]
for index, row in from_cd_merged.iterrows():
    if str(row["direction"])==str(1):
        strand.append("+")
    if str(row["direction"])==str(-1):
        strand.append("-")
    else:
        pass
from_cd_merged["strand"] = strand
from_cd_merged=from_cd_merged.rename(columns={"first_x" : "first", "last_x" : "last"})
# from_cd_merged["first"] = first
# from_cd_merged["last"] = last
from_cd_merged["frame"] = 0
dummy=[]
dummy2=[]
for index,row in from_cd_merged.iterrows():
    dummy.append( "gene_id \"KH" + str(index) + "\"; transcript_id \"KH" + str(index) + ".1\"")
    #dummy2.append(str(row["names"])+"_"+str(row["chr"]) + ":" + str(row["first"]) + "-" + str(row["last"]) + "-" + str(int(float(row["direction"]))))
    dummy2.append(str(row["names"])+"_"+str(row["chr"]) + ":" + str(int(float(str(row["first"]).replace(".0","").replace("--1", "").replace("-1", "")))) + "-" + str(int(float(str(row["last"]).replace(".0","").replace("--1", "").replace("-1", "")))) + "-" + str(int(float(row["direction"]))))
from_cd_merged["attribute"] = dummy
from_cd_merged["query"] = dummy2


from_cd_merged.to_csv(BLAST_species+"_cdhit_full.csv")

from_cd_merged_rearranged = from_cd_merged[["chr", "source", "feature", "first", "last", "score", "strand", "frame", "attribute"]]

from_cd_merged_rearranged.to_csv (BLAST_species + "_from_pipeline.gtf", sep="\t", header=False, index=False,quoting=csv.QUOTE_NONE)

from_cd_forrecip = from_cd_merged [["names", "chr", "first", "last", "seqs", "direction", "frame", "strand","attribute","query"]]

pickle_out=open("from_cd_forrecip.pickle","wb")
pickle.dump(from_cd_forrecip, pickle_out)
pickle_out.close()
