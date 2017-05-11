#!/usr/bin/env Rscript

suppressMessages(library(gdata)) # contains write.fwf
source("scripts/aux/functions.r")

args = commandArgs(trailingOnly=TRUE)

stopifnot(length(args) == 5)

# input files
file_loci_desc  = args[1]
file_taxa       = args[2]
file_filter_cfg = args[3]

# output files
file_taxa_desc   = args[4]
dir_filter_stats = args[5]

# output graphs
img_taxa_per_locus=paste(dir_filter_stats, "taxa_per_locus.eps", sep='/')
img_taxa_per_locus_filtered=paste(dir_filter_stats, "taxa_per_locus_filtered.eps", sep='/')
img_varsites_per_locus=paste(dir_filter_stats, "varsites_per_locus.eps", sep='/')
img_varsites_per_locus_filtered=paste(dir_filter_stats, "varsites_per_locus_filtered.eps", sep='/')
img_gaps_per_locus=paste(dir_filter_stats, "gaps_per_locus.eps", sep='/')
img_gaps_per_locus_filtered=paste(dir_filter_stats, "gaps_per_locus_filtered.eps", sep='/')
img_taxa_rep=paste(dir_filter_stats, "img_taxa_representation.eps", sep='/')

# read taxa names
tnames = scan(file_taxa, "", quiet=TRUE)

# read loci description
data   = parse_loci_desc(file_loci_desc)

ntaxa  = length(tnames)
nloci  = length(data$id)

# convert taxa map into a matrix
taxamap = matrix(
            unlist(
              lapply(
                data$tmap, 
                function(x) as.numeric(strsplit(x, split="")[[1]]))), 
            ncol=ntaxa, 
            byrow=TRUE)

taxadata = list()
taxadata$names = tnames
taxadata$nloci   = sapply(1:ntaxa, function(x) sum(taxamap[,x]))
taxadata$percent = 100 * taxadata$nloci / nloci

#filter criteria
filter = parse_filter(file_filter_cfg, data)

nloci.filtered = sum(filter$all)
taxadata$nloci.filtered   = sapply(1:ntaxa, function(x) sum(taxamap[which(filter$all==1),x]))
taxadata$percent.filtered = 100 * taxadata$nloci.filtered / nloci.filtered

cat("Filtered data contains", sum(filter$all), "loci\n")

taxamatrix = matrix(unlist(taxadata), ncol=length(taxadata), byrow=FALSE)
colnames(taxamatrix) = names(taxadata)

write(x = paste("#", paste(names(taxadata), collapse=" ")), file = file_taxa_desc, append = FALSE)
write.fwf(
  x     = taxamatrix,
  file  = file_taxa_desc,
  quote = FALSE,
  rownames = FALSE,
  colnames = FALSE,
  append = TRUE)
cat(">> Taxa description written to", file_taxa_desc, "\n")

# plot graphs

setEPS()

# colors = c("red"), "yellow", "green", "violet", "orange", "blue", "pink", "cyan")

postscript(img_taxa_per_locus)
hist(100*data$tprop, breaks=25, col=c("cyan"), xlab="% of taxa", main="Percentage of taxa per locus")
invisible(dev.off())
postscript(img_taxa_per_locus_filtered)
hist(100*data$tprop[filter$all==1], breaks=25, col=c("cyan"), xlab="% of taxa", main="Percentage of taxa per locus")
invisible(dev.off())
cat(">> Histogram taxa/locus dumped to", img_taxa_per_locus, "\n")

postscript(img_varsites_per_locus)
hist(100*data$vprop, breaks=25, col=c("cyan"), xlab="% of variable sites", main="Percentage of variable sites per locus")
invisible(dev.off())
postscript(img_varsites_per_locus_filtered)
hist(100*data$vprop[filter$all==1], breaks=25, col=c("cyan"), xlab="% of variable sites", main="Percentage of variable sites per locus")
invisible(dev.off())
cat(">> Histogram variable sites/locus dumped to", img_varsites_per_locus, "\n")

postscript(img_gaps_per_locus)
hist(100*data$gapy, breaks=25, col=c("cyan"), xlab="% of gaps", main="Gapyness per locus")
invisible(dev.off())
postscript(img_gaps_per_locus_filtered)
hist(100*data$gapy[filter$all==1], breaks=25, col=c("cyan"), xlab="% of gaps", main="Gapyness per locus")
invisible(dev.off())
cat(">> Histogram gaps/locus dumped to", img_gaps_per_locus, "\n")

cat("\n  %taxa    n_loci   %loci    n_loci[F]   %loci[F]\n")
cat("--------------------------------------------------\n")

for (i in 0:9)
{
  n_samples = sum(data$tprop >= i/10)
  n_samples_filtered = sum(data$tprop[filter$all==1] >= i/10)
  cat(">=", formatC(i*10, width=5), "% : ", 
      formatC(n_samples, width=6), " : ", 
      format(100*n_samples/nloci, nsmall=3, digits=3, width=8, justify="right"), "%", 
      formatC(n_samples_filtered, width=8),
      format(100*n_samples_filtered/sum(filter$all), nsmall=3, digits=3, width=8, justify="right"), "%\n", 
      sep="")
}
cat("\n")

cat("\n  ftaxa    n_loci   %loci    n_loci[F]   %loci[F]\n")
cat("--------------------------------------------------\n")
for (i in sort(unique(data$eftaxa)))
{
  n_samples = sum(data$eftaxa >= i)
  n_samples_filtered = sum(data$eftaxa[filter$all==1] >= i)
  cat(">=", formatC(i, width=4), " : ", 
      formatC(n_samples, width=6), " : ", 
      format(100*n_samples/nloci, nsmall=3, digits=3, width=8, justify="right"), "%", 
      formatC(n_samples_filtered, width=8),
      format(100*n_samples_filtered/sum(filter$all), nsmall=3, digits=3, width=8, justify="right"), "%\n", 
      sep="")
}
cat("\n")

taxamap = lapply(strsplit(data$tmap, ''), as.numeric)
n_taxa = length(taxamap[[1]])

postscript(img_taxa_rep)
hist(taxadata$percent, breaks=25, col=c("orange"), xlab="Coverage (%)", main="Taxa representation")
invisible(dev.off())
cat(">> Histogram taxa representation dumped to", img_taxa_rep, "\n")
