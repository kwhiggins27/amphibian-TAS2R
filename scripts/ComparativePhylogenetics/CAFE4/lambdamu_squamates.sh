#!shell
date
version

#specify data file, p-value threshold, # of threads to use, and log file
load -i squamates_tas2r_numbers.fixedlabels.txt -t 12 -l squamates.txt

#the phylogenetic tree structure with branch lengths
tree ((Euleptes-europaea:111.71,Eublepharis-macularius:111.71):77.58994,(Hemicordylus-capensis:173.55365,(((Ahaetulla-prasina:53.33742,((((Bungarus-multicinctus:29.14278,Hydrophis-cyanocinctus:29.14278):2.65239,Naja-naja:31.79517):6.06175,Thamnophis-elegans:37.85692):0.34278,((Vipera-ursinii:12.78732,Vipera-latastei:12.78732):20.84418,Crotalus-viridis:33.6315):4.5682):15.13772):107.70799,(Anolis-carolinensis:84.83432,(Phrynosoma-platyrhinos:54.11697,(Sceloporus-undulatus:28.76711,Sceloporus-tristichus:28.76711):25.34986):30.71735):76.21109):5.85832,((((Podarcis-raffonei:17.60364,Podarcis-muralis:17.60364):25.69284,Lacerta-agilis:43.29648):3.60488,Zootoca-vivipara:46.90136):79.85864,Rhineura-floridana:126.76):40.14373):6.64992):15.74629)

#search for 2 parameter model
lambdamu -s


