library(ape)
library(phytools)
library(TreeTools)


## Read in tree

#tr=read.tree("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Oct_Tree3_5amp_nt/5core_nucleic_acid_fish.fas.treefile")
tr=read.tree("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA/combo_aa.aln.treefile")
#tr=read.tree("/Volumes/wengpj01/phylogenetics/renamed_and_numbered8.tree")


annotaiton_raw="/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Nov_aminoMflag/combo_Mflag.csv"
#annotation_raw="/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Nov_amino/combo_noflag.csv"
annotation=read.csv(annotation_raw)
#create vector for colors with the same number of elements as tips on the tree

cols=rep("black",length(tr$tip.label))

#Now use grep to find which tips match each group and set the corresponding elements in the color vector
#When we need to grep for more than one pattern (e.g. frogs are Amphibia, Xetr, Rhma, etc....) we can specify them as "P1|P2|P3"

##NEED TO COMPLETE CATEGORIES HERE"

# amphibians=grep("Amphibia|Xetr|Rhma|Amme", tr$tip.label)
# cols[amphibians]="darkgreen"

axolotl = grep("axolotl", tr$tip.label)
cols[axolotl]="#D55E00" #"darkgreen"

bullfrog = grep("bullfrog", tr$tip.label)
cols[bullfrog]="#166A53"#"olivedrab"

cane = grep("cane", tr$tip.label)
cols[cane]="#œ" #lightgreen"

dart=grep("dart", tr$tip.label)
cols[dart]="#063781" #"red"

clawed=grep("clawed", tr$tip.label)
cols[clawed]="#AEF2B5" #CBA9F4"

outgroup=grep("ora",tr$tip.label)
cols[outgroup]="#9A9F9B" #"plum4"

turt=grep("Testudines", tr$tip.label)
cols[turt]="#9245E3" #"purple1"


# #pdf(file="TipTree_TASR.pdf", w=12, h=12 , useDingbats=F)
# pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Oct_Tree3_5amphibian/5core_amino_acid.pdf")
# plot(tr, type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
# dev.off()
# #dev.off()

pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_numbered_6.pdf")
plot(tr, type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
ape::nodelabels()
dev.off()

tr5 <- reroot(tr, 644)
ape::write.tree(tr5, file="/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/reroot.newick")

pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_numbered_7.pdf")
plot(tr5, type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
ape::nodelabels()
dev.off()

is_tip <- tr5$edge[,2] <= length(tr5$tip.label)
ordered_tips <- rev(tr5$tip.label)




# Create an empty matrix called reordered
reordered <- data.frame(matrix(NA, nrow = length(ordered_tips), ncol = ncol(annotation)))
# 
# # Set the first column of reordered as ordered list
# reordered[, 1] <- unlist(ordered)

# Set the column names of reordered
colnames(reordered) <- colnames(annotation)

# Loop through each value in ordered
for (i in 1:length(ordered_tips)) {
  # Search for the value in ordered within data
  match_row <- annotation[annotation$X == ordered_tips[i], ]
  
  # Check if a match was found
  if (nrow(match_row) == 1) {
    # Append the matched row to reordered
    reordered[i, ] <- match_row
  } else {
    # Handle the case when no match is found
    # You can choose to do nothing or perform some other action here
    # For this example, I'll set all columns to NA for the current row
    reordered[i, -1] <- NA
  }
}



reordered[, 2:8] <- lapply(reordered[, 2:8], function(x) ifelse(is.na(x), 0, as.numeric(x)))

matrix_data=reordered[,2:8]
# Define a small positive value to add to the data
epsilon <- 1e-1  # You can adjust the value as needed

# Take the logarithm of each value in matrix_data
log_matrix_data <- log(matrix_data + epsilon)




rownames(log_matrix_data) <- ordered_tips
pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/heatmap_aa.pdf")
phylo.heatmap(tr,log_matrix_data,standardize=FALSE)
dev.off()

node_labels <-tr5$node.label
node_labels <- gsub("^/", "", node_labels)  # Remove leading fowardslash
node_labels <- as.numeric(node_labels)  # Convert the cleaned strings to numeric
node_labels[is.na(node_labels)] <- 1

# Custom function to remove leading '/' and convert to numeric
clean_and_convert <- function(x) {
  x <- gsub("^/", "", x)  # Remove leading '/'
  as.numeric(x)           # Convert to numeric
}
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

tr5 <- reroot(tr, 644)
pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_2.pdf")
plot(tr5, type="phylogram", TRUE, cex=0.15, tip.color=cols, lwd=0.3)
ape::nodelabels()
dev.off()



#pdf(file="TipTree_TASR.pdf", w=12, h=12 , useDingbats=F)
pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_2.pdf")
plot(tr5, type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
nodelabels(pch = 21, col = "black", cex = 0.4, frame = "n", bg = sapply(node_labels, fill_color))
dev.off()
#dev.off()

pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_no_nodes_2.pdf")
plot(tr5, type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
dev.off()



# tr2=Preorder(tr)
# # tr2=reorder(tr, order = "cladewise", index.only = FALSE)
# # cladewise(tr2)
# # postorder(tr2)
# 
# pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_numbered_6.pdf")
# plot(tr2, type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
# ape::nodelabels()
# dev.off()
# #
#
# ##Figure 5E
# pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_numbered_8.pdf")
# plot(Subtree(Subtree(tr2,801),130), type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
# ape::nodelabels()
# dev.off()
# 
# tr3=Subtree(Subtree(tr2,801),130)
# ape::write.tree(tr3, file="/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5E_subset.newick")
# 
# ##Figure 5D
# pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_numbered_8.pdf")
# plot(Subtree(Subtree(tr2,1038),150), type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
# ape::nodelabels()
# dev.off()
# 
# tr3=Subtree(Subtree(tr2,1038),150)
# ape::write.tree(tr3, file="/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5D_subset.newick")
# 
# ##Figure 5C
# pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_numbered_8.pdf")
# plot(Subtree(Subtree(tr2,740),40), type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
# ape::nodelabels()
# dev.off()
# 
# tr3=Subtree(Subtree(tr2,740),40)
# ape::write.tree(tr3, file="/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5C_subset.newick")
# 
# 
# #Figure 5B
# pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_numbered_8.pdf")
# plot(Subtree(Subtree(tr2,793), 77), type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
# ape::nodelabels()
# dev.off()
# 
# tr3=Subtree(Subtree(tr2,793), 77)
# ape::write.tree(tr3, file="/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5B_subset.newick")

# pdf("/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Jan_5amp_AA_MM/5core_amino_acid_nodes_numbered_8.pdf")
# plot(Subtree(Subtree(tr2, 676),276), type="phylogram", FALSE, cex=0.15, tip.color=cols, lwd=0.3)
# ape::nodelabels()
# dev.off()

# tr3=Subtree(Subtree(tr2, 676),276)
# ape::write.tree(tr3, file="/Users/kwhiggins27/Desktop/OneDrive - Harvard University/Weng Lab/TAS2R evolution/Nov_amino/chicken_subset.newick")

