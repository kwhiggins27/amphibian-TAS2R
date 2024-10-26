#PGLS on multiple traits related to number of TAS2R

library(ape)
library(geiger)
library(phylolm)

tree=read.tree("Timetree_ultrametric.tre")
data_trim=read.delim("tas2r_Github_data.tsv", stringsAsFactors=F)

tree=drop.tip(tree, name.check(tree,data_trim)$tree_not_data)

name.check(tree,data_trim)


# Phylogenetic gls

#Number of genes vs genome size
summary(phylolm(log(Number.of.Genes.x+1)~log(genome_size), data=data_trim, tree,model = "OUfixedRoot")) # beta=0.2308, t=1.60, p=0.11)

## Now just amphibians
summary(phylolm(log(Number.of.Genes.x+1)~log(genome_size), data=data_trim[data_trim$PlottingClade%in%c("Caudata","Anura","Gymnophiona"),], tree,model = "OUfixedRoot")) #beta=0.132, t=0.429, p=0.67


## Now clusters and gene counts

## Full model
summary(phylolm(log(Number.of.Genes.x+1)~log(clusters+1)+log(genes_per_cluster+1), data=data_trim, tree, model = "OUfixedRoot"))

## Just batrachians
summary(phylolm(log(Number.of.Genes.x+1)~log(clusters+1)+log(genes_per_cluster+1), data=data_trim[data_trim$PlottingClade%in%c("Caudata","Anura"),], tree, model = "OUfixedRoot"))

## Last thing: Model including genome size for the supplement

summary(phylolm(log(Number.of.Genes.x+1)~log(clusters+1)+log(genes_per_cluster+1)+log(genome_size), data=data_trim, tree, model = "OUfixedRoot"))
