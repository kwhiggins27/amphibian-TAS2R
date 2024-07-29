library(ape)

## Read in tree

setwd("/Volumes/wengpj01/amphibian-TAS2R/results/phylogenetics") #Modify this to your actual working directory
tr=read.tree("../../results/phylogeneticsNovtree_skeleton_10K.fas.treefile")  #Modify.  This tree will need to be separately download from the manuscript's webpage

cols=rep("black",length(tr$tip.label))

accession_1="GCA_009819655.1" #Modify.  One of the test organisms
accession_2="GCA_014839755.1" #Modify.  One of the test organisms
accession_3="GCA_015227805.3" #Modify.  Outgroup

out_file_nodes="crows.pdf" #Modify
out_file_no_nodes="crows_no_nodes.pdf" #Modify

accessions <- c(accession_1, accession_2, accession_3) #Put the three desired accession numbe4rs here and in lines 16, 20, 23
# Subset tip labels based on patterns
subset_labels <- tr$tip.label[startsWith(tr$tip.label, accessions[1]) | startsWith(tr$tip.label, accessions[2]) | startsWith(tr$tip.label, accessions[3])]


tr2 = keep.tip(tr, subset_labels)
common = grep(accession_1, tr2$tip.label) 
cols[common]="#D55E00"

wood = grep(accession_2, tr2$tip.label)
cols[wood]="#166A53"

pyx = grep("accession_3", tr2$tip.label)
cols[pyx]="#063781" 


# Define a function to determine the fill color of the circles based on node.label values
fill_color <- function(x) {
  if (x > 0.95) {
    return("red")  # Greater than 0.95, make the fill red
  } else if (x >= 0.90) {
    return("orange")         # 0.9 to 0.95, make the fill orange
  } else {
    return("yellow")        # Less than 0.9, make the fill yellow
  }
}

node_labels <-tr2$node.label
node_labels <- gsub("^/", "", node_labels)  # Remove leading fowardslash
node_labels <- as.numeric(node_labels)  # Convert the cleaned strings to numeric
node_labels[is.na(node_labels)] <- 1

pdf(out_file_nodes)
plot(tr2, type = "phylogram", FALSE, cex = 0.15,lwd = 0.3,tip.color=cols)
nodelabels(pch = 21, cex = 0.5,frame = "n", bg = sapply(node_labels, fill_color))
dev.off()

pdf(out_file_no_nodes)
plot(tr2, type = "phylogram", FALSE, cex = 0.15,lwd = 0.3,tip.color=cols)
dev.off()

