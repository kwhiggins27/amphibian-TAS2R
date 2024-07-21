library(ape)

## Read in tree

setwd("/Volumes/wengpj01/phylogenetics/orth_tree")
tr=read.tree("/Volumes/wengpj01/phylogenetics/Nov_bigtree/Novtree_skeleton_10K.fas.treefile")

cols=rep("black",length(tr$tip.label))

accessions <- c("GCA_009819655.1", "GCA_014839755.1", "GCA_015227805.3")
# Subset tip labels based on patterns
subset_labels <- tr$tip.label[startsWith(tr$tip.label, accessions[1]) | startsWith(tr$tip.label, accessions[2]) | startsWith(tr$tip.label, accessions[3])]


tr2 = keep.tip(tr, subset_labels)
common = grep("GCA_009819655.1", tr2$tip.label)
cols[common]="#D55E00"

wood = grep("GCA_014839755.1", tr2$tip.label)
cols[wood]="#166A53"

pyx = grep("GCA_015227805.3", tr2$tip.label)
cols[pyx]="#063781" 


# Define a function to determine the fill color of the circles based on node.label values
fill_color <- function(x) {
  if (x > 0.95) {
    return("red")  # Greater than 0.5, make the circle transparent
  } else if (x >= 0.90) {
    return("orange")         # 0.25 to 0.5, make the fill gray
  } else {
    return("yellow")        # Less than 0.25, make the fill black
  }
}

node_labels <-tr2$node.label
node_labels <- gsub("^/", "", node_labels)  # Remove leading fowardslash
node_labels <- as.numeric(node_labels)  # Convert the cleaned strings to numeric
node_labels[is.na(node_labels)] <- 1

pdf("crows.pdf")
plot(tr2, type = "phylogram", FALSE, cex = 0.15,lwd = 0.3,tip.color=cols)
nodelabels(pch = 21, cex = 0.5,frame = "n", bg = sapply(node_labels, fill_color))
dev.off()

pdf("crows_no_nodes.pdf")
plot(tr2, type = "phylogram", FALSE, cex = 0.15,lwd = 0.3,tip.color=cols)
dev.off()

