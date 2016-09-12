descfile="loci.desc"
taxafile="taxa.desc"

data=scan(descfile, list(id=0, ntax=0, tprop=0, len=0, nvar=0, ninf=0, vprop=0, gapy=0, tmap=""), comment.char="#")
taxadata=scan(taxafile, list(name="", nloci=0, percent=0), comment.char="#")

hist(100*data$tprop, breaks=25, col=c("cyan"), xlab="% of taxa", main="Percentage of taxa per locus")
hist(100*data$vprop, breaks=25, col=c("cyan"), xlab="% of variable sites", main="Percentage of variable sites per locus")
hist(100*data$gapy, breaks=25, col=c("cyan"), xlab="% of gaps", main="Gapyness per locus")

n_loci = length(data$id)

for (i in 0:9)
{
  n_samples = sum(data$tprop >= i/10)
  cat(">=", i/10, "% : ", n_samples, " : ", n_samples/n_loci, "%\n", sep="")
}

taxamap = lapply(strsplit(data$tmap, ''), as.numeric)
n_taxa = length(taxamap[[1]])

filter = which(data$tprop > 0.74)
msataxa = numeric(n_taxa)
for (i in filter)
  msataxa = msataxa + taxamap[[i]]
comp_taxa = sum(msataxa == length(filter))
cat ("Comprehensive taxa = ", comp_taxa, "\n")

for (j in 1:n_taxa)
{
  msataxon = numeric(n_taxa)
  for (i in filter)
    if (taxamap[[i]][j] == 1)
      msataxon = msataxon + taxamap[[i]]
  if (sum(msataxon == 0) > 0)
    cat (j, "  :  ", which(msataxon == 0), "\n")
}
# colors = c("red"), "yellow", "green", "violet", "orange", "blue", "pink", "cyan")

hist(taxadata$percent, breaks=25, col=c("orange"), xlab="Coverage (%)", main="Taxa representation")
