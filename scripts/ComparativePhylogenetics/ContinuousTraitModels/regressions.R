
## l1OU and PGLS of TAS2R gene number

library(ape)
library(geiger)
library(phytools)
library(ips)
library(l1ou)
library(phylolm)
library(mjcbase)

#Read in stuff
tree=read.tree("Timetree_ultrametric.tre")
gen_data=read.csv("genes_clusters_etc_removedFailed.csv", stringsAsFactors=F)
rownames(gen_data)=gen_data$latin.x

## QUickly transform the number of clusters from NA to 0 if there are more than 2 genes

for(i in 1:nrow(gen_data)){
	if(is.na(gen_data$clusters[i]) && gen_data$Number.of.Genes.x[i]>1){gen_data$clusters[i]=0}
}

genomeSize=read.csv("all_genome_sizes.csv", stringsAsFactors=F)
rownames(genomeSize)=genomeSize$latin

## Merge both tables (and get just the variables we want)

data=cbind.smart(gen_data[,c(1,4,7,8:10,13,23)],genomeSize[,5:6])

## Namecheck

remove=which(rownames(data)%in%name.check(tree,data)$data_not_tree)
data_trim=data[-remove,]
tree=drop.tip(tree, name.check(tree,data)$tree_not_data)

name.check(tree,data_trim)

summary(lm(log(Number.of.Genes.x+1)~log(genome_size), data=data_trim)) # beta=1.08, t=12.25, p<2e-16

# Phylogenetic gls
summary(phylolm(log(Number.of.Genes.x+1)~log(genome_size), data=data_trim, tree,model = "OUfixedRoot")) # beta=0.2308, t=1.60, p=0.11)

## Now just amphibians
summary(phylolm(log(Number.of.Genes.x+1)~log(genome_size), data=data_trim[data_trim$PlottingClade%in%c("Caudata","Anura","Gymnophiona"),], tree,model = "OUfixedRoot")) #beta=0.132, t=0.429, p=0.67

## Standalone plot for this
par(mfrow=c(1,2))
plot(data_trim$genome_size, data_trim$Number.of.Genes, pch=21, bg=colors, xlab="Genome Size (c-value)", ylab="Number of Genes", cex.lab=1.1, cex=1.5, main="All Vertebrates")
#legend(25, 250, c(""))

## Just amphibians
plot(data_trim[data_trim$PlottingClade%in%c("Caudata","Anura","Gymnophiona"),]$genome_size, data_trim[data_trim$PlottingClade%in%c("Caudata","Anura","Gymnophiona"),]$Number.of.Genes, pch=21, bg=colors[data_trim$PlottingClade%in%c("Caudata","Anura","Gymnophiona")], xlab="Genome Size (c-value)", ylab="Number of Genes", cex.lab=1.1, cex=1.5, main="Amphibians only")

## Now clusters and gene counts
# First get counts instead of fractions for a couple values (NOT USED, these correlations don't make a lot of senseto do):

#data_trim$NumberClustered=round(data_trim$Fraction.Clustered*data_trim$Number.of.Genes.x)
#data_trim$NumberNearest_in_cluster=round(data_trim$Nearest_in_cluster*data_trim$Number.of.Genes.x)

## Full model
summary(phylolm(log(Number.of.Genes.x+1)~log(clusters+1)+log(genes_per_cluster+1), data=data_trim, tree, model = "OUfixedRoot"))

## Just batrachians
summary(phylolm(log(Number.of.Genes.x+1)~log(clusters+1)+log(genes_per_cluster+1), data=data_trim[data_trim$PlottingClade%in%c("Caudata","Anura"),], tree, model = "OUfixedRoot"))

summary(lm(log(Number.of.Genes.x+1)~log(genome_size), data=data_trim)) # beta=1.08, t=12.25, p<2e-16

# Phylogenetic gls
summary(phylolm(log(Number.of.Genes.x+1)~log(genome_size), data=data_trim, tree,model = "OUfixedRoot")) # beta=0.2308, t=1.60, p=0.11)

## Now just amphibians
summary(phylolm(log(Number.of.Genes.x+1)~log(genome_size), data=data_trim[data_trim$PlottingClade%in%c("Caudata","Anura","Gymnophiona"),], tree,model = "OUfixedRoot")) #beta=0.132, t=0.429, p=0.67

## Standalone plot for this
par(mfrow=c(1,2))
plot(data_trim$genome_size, data_trim$Number.of.Genes, pch=21, bg=colors, xlab="Genome Size (c-value)", ylab="Number of Genes", cex.lab=1.1, cex=1.5, main="All Vertebrates")
#legend(25, 250, c(""))

## Just amphibians
plot(data_trim[data_trim$PlottingClade%in%c("Caudata","Anura","Gymnophiona"),]$genome_size, data_trim[data_trim$PlottingClade%in%c("Caudata","Anura","Gymnophiona"),]$Number.of.Genes, pch=21, bg=colors[data_trim$PlottingClade%in%c("Caudata","Anura","Gymnophiona")], xlab="Genome Size (c-value)", ylab="Number of Genes", cex.lab=1.1, cex=1.5, main="Amphibians only")

## Now clusters and gene counts
# First get counts instead of fractions for a couple values (NOT USED, these correlations don't make a lot of senseto do):

#data_trim$NumberClustered=round(data_trim$Fraction.Clustered*data_trim$Number.of.Genes.x)
#data_trim$NumberNearest_in_cluster=round(data_trim$Nearest_in_cluster*data_trim$Number.of.Genes.x)

## Full model
summary(phylolm(log(Number.of.Genes.x+1)~log(clusters+1)+log(genes_per_cluster+1), data=data_trim, tree, model = "OUfixedRoot"))

## Just batrachians
summary(phylolm(log(Number.of.Genes.x+1)~log(clusters+1)+log(genes_per_cluster+1), data=data_trim[data_trim$PlottingClade%in%c("Caudata","Anura"),], tree, model = "OUfixedRoot"))


