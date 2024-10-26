#!shell
date
version

#specify data file, p-value threshold, # of threads to use, and log file
load -i batrachians_tas2r_numbers.fixedlabels.txt -t 12 -l batrachians.txt

#the phylogenetic tree structure with branch lengths
tree ((Ambystoma-mexicanum:156.21,Pleurodeles-waltl:156.21):135.21676,((((Spea-bombifrons:137.8525,(Leptobrachium-leishanense:9.5992,Leptobrachium-ailaonicum:9.5992):128.2533):53.8805,((((((Bufo-bufo:13.2894,Bufo-gargarizans:13.2894):54.74655,Engystomops-pustulosus:68.03595):4.57428,(Hyla-sarda:54.3227,Dendropsophus-ebraccatus:54.3227):18.28753):13.87485,Eleutherodactylus-coqui:86.48508):39.90536,Pseudophryne-corroboree:126.39044):22.71414,((Pyxicephalus-adspersus:82.9102,(Lithobates-sylvaticus:33.43947,((Rana-kukunoris:19.24698,Rana-temporaria:19.24698):0.85122,Rana-muscosa:20.0982):13.34127):49.47073):29.0498,Gastrophryne-carolinensis:111.96):37.14458):42.62842):10.60508,(Hymenochirus-boettgeri:108.08258,((Xenopus-borealis:39.41876,Xenopus-laevis:39.41876):18.81114,Xenopus-tropicalis:58.2299):49.85268):94.2555):7.37192,Bombina-bombina:209.71):81.71676)

#search for one parameter model
#lambda -s

#search for 2 parameter model
lambdamu -s


