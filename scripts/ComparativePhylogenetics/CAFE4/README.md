## Running CAFE4 on individual vertebrate lineages

CAFE4 was run 50 times on each input file using a simple loop, which saved the last line of the output file, containing the model parameters with the highest likelihood for each run, to a text file. Below is an example using the squamate data.

``
for i in {1..50}; do cafe lambdamu_squamates.sh; mv squamates.txt squamates_"$i".txt; tail -n 1 squamates_"$i".txt >> squamates_runs.txt; done
``
<br>
<br>
We then extracted the run with the highest likelihood in R:

```R
files=c("birds_runs.txt","squamates_runs.txt", "batrachians_runs.txt", "RFfish_runs.txt")

lambda=c()
mu=c()
lnL=c()

for(i in 1:length(files)){
	tab=read.table(files[i])
	bestRun=tab[which.min(tab$V11),]
	lambda[i]=bestRun$V3[1]
	mu[i]=bestRun$V8[1]
	lnL[i]=-1*bestRun$V11[1]
	}
res=data.frame(files, lambda, mu, lnL)
write.table(res, "cafe4_rates.txt")
```
The output file with the parameter estimates for each clade is also in the repository.
