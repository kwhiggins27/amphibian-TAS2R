#PGLS on multiple traits related to number of TAS2R

library(ape)
library(geiger)
library(phylolm)
library(l1ou)
library(phytools)

# Read in data

tree=read.tree("Timetree_ultrametric.tre")
data_trim=read.delim("tas2r_Github_data.tsv", stringsAsFactors=F)

# Rotate nodes for plotting
tree=rotateNodes(tree,c(664,939,966))

## Name-check
tree=drop.tip(tree, name.check(tree,data_trim)$tree_not_data)
name.check(tree,data_trim)

#Create a variable with log-transformed family size
log_num=log(data_trim$Number.of.Genes+1)
names(log_num)=rownames(data_trim)


# Run l1ou

data_adj=adjust_data(tree, log_num, normalize=F)

#fullFit=estimate_shift_configuration(data_adj$tree, data_adj$Y, max.nShifts=50, criterion="pBIC", nCores = 16) #Multicore doesn't report back models for some reason

fullFit=estimate_shift_configuration(data_adj$tree, data_adj$Y, max.nShifts=50, criterion="pBIC")

## Make table with pBIC weights

l1ou_mods=data.frame(Branches=as.character(fullFit$profile$configurations), pBIC=fullFit$profile$scores)

for (i in 1:nrow(l1ou_mods)){
	l1ou_mods$d_pBIC[i]=l1ou_mods$pBIC[i]-min(l1ou_mods$pBIC)
	l1ou_mods$wt_pBIC[i]=exp(-0.5*l1ou_mods$d_pBIC[i])
	}
l1ou_mods$wt_pBIC=l1ou_mods$wt_pBIC/sum(l1ou_mods$wt_pBIC)

# check out the top five models

head(l1ou_mods, 5)


## Now do more focused model fitting

## Best two l1ou models excluding Setophaga coronata for simplicity. Note we're scoring based on AIC now


OU_6shift_nb=fit_OU(fullFit$tree, fullFit$Y, c(257, 677, 688, 740, 1270, 1282), criterion="AIC")
OU_5shift=fit_OU(fullFit$tree, fullFit$Y, c(257, 677, 740, 1270, 1282), criterion="AIC")

## Now create simmap objects with the different regimes for brownie.lite. This needs to be done branch by branch

fiveshift=paintSubTree(tree, node=655, state="1", anc="0")
fiveshift=paintSubTree(fiveshift, node=918, state="2")
fiveshift=paintSubTree(fiveshift, node=923, state="3")
fiveshift=paintSubTree(fiveshift, node=1162, state="4")
fiveshift=paintSubTree(fiveshift, node=952, state="5")

sixshift_nb=paintSubTree(fiveshift, node=966, state="6")

## Make sure they look good

plot(fiveshift)
plot(sixshift)

# Now fit brownian models. Increasing maxit was necessary to reach convergence

BM_5shift=brownie.lite(fiveshift, log_num, maxit=10000)
BM_6shift_nb=brownie.lite(sixshift_nb, log_num, maxit=10000)

## Finally fit single-regime BM and OU

BM=fitContinuous(tree, log_num, model="BM", control=list(niter=500))
OU=fitContinuous(tree, log_num, model="OU", control=list(niter=500))
lambda=fitContinuous(tree, log_num, model="lambda", control=list(niter=500))

## Make table with AIC scores and weights

modComp_nb=data.frame(Model=c("OU 6-shift", "OU 5-shift", "BM 6-shift", "BM 5-shift", "OU_single", "BM_single","Lambda"), AIC=c(OU_6shift_nb$score, OU_5shift$score, AIC(BM_6shift_nb)[2], AIC(BM_5shift)[2], OU$opt$aic, BM$opt$aic,lambda$opt$aic))

modComp_nb=modComp_nb[order(modComp_nb$AIC),]


for (i in 1:nrow(modComp_nb)){
	modComp_nb$dAIC[i]=modComp_nb$AIC[i]-min(modComp_nb$AIC)
	modComp_nb$wt_AIC[i]=exp(-0.5*modComp_nb$dAIC[i])
	}
	
modComp_nb$wt_AIC=modComp_nb$wt_AIC/sum(modComp_nb$wt_AIC)

# Done
