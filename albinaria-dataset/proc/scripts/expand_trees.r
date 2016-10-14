#!/usr/bin/env Rscript

# Usage:
# Rscript --vanilla 

args = commandArgs(trailingOnly=TRUE)

#input files
taxafile = "taxa"
taxamapfile = "taxamap"
treefile = args[1]
stopifnot(file.exists(treefile))
outputfile = paste(treefile, "expanded", sep=".") #"consensus.trees.expanded"
locifile = "loci.desc"

stopifnot(file.exists(locifile))
stopifnot(file.exists(taxafile))

if (file.exists(taxamapfile)) {
  taxamap = read.table(file=taxamapfile)
} else {
  data = scan(locifile, list(id=0, ntax=0, tprop=0, len=0, nvar=0, ninf=0, vprop=0, gapy=0, tmap="", eftaxa=0, dups=""), comment.char="#", quiet=TRUE)

  taxamap = lapply(strsplit(data$tmap, ''), as.numeric)
  n_taxa = length(taxamap[[1]])

  filterInf = which(data$ninf == 0)
  filterTax = which(data$eftaxa < 4)
  filterOut = unique(c(filterInf, filterTax))
  
  data_filtered = lapply(data, function(x) x[-filterOut])
  taxamap_filtered = taxamap[-filterOut]

  # update duplicates
  seqduplicates = lapply(strsplit(data_filtered$dups, ','), as.numeric)
  stopifnot(length(seqduplicates) == length(taxamap_filtered))
  for (i in 1:length(seqduplicates))
  {
    duplist = seqduplicates[[i]]
    for (j in 0:((length(duplist)/2)-1))
    {
      if (length(duplist)%%2 == 0)
      {
        i1 = duplist[j*2+1]
        i2 = duplist[j*2+2]
        if (i1 != i2)
        {
          x1 = taxamap_filtered[[i]][i1]
          x2 = taxamap_filtered[[i]][i2]
          nextV = max(taxamap_filtered[[i]]) + 1
          if (x1 == 1 && x2 == 1)
          {
            taxamap_filtered[[i]][i1] = taxamap_filtered[[i]][i2] = nextV
          }
          else
          {
            if (x1 > 1)
            {
              taxamap_filtered[[i]][i2] = taxamap_filtered[[i]][i1]
            }
            else
            {
              taxamap_filtered[[i]][i1] = taxamap_filtered[[i]][i2]
            }
          }
        }
      }
    }
  }
  tmpeftaxa = sapply(taxamap_filtered, function(x) sum(x==1) + length(unique(x)) - 2)

  taxamap_table=matrix(nrow=length(taxamap_filtered),
                       ncol=length(taxamap_filtered[[1]]) + 1)
  for (i in 1:length(taxamap_filtered))
  {
    taxamap_table[i,1]  = data_filtered$id[i]
    taxamap_table[i,2:(n_taxa+1)] = taxamap_filtered[[i]]
  }

  taxamap = taxamap_table
}

trees = scan(treefile, "", quiet=TRUE)
taxanames = scan(taxafile, "", quiet=TRUE)
n_trees = length(trees)
row_len = dim(taxamap)[2]

stopifnot(dim(taxamap)[1] == n_trees)
stopifnot(length(taxanames) == row_len - 1)

perc_step = n_trees %/% 100
perc_iter = perc_step
perc_count = 0
for (i in 1:n_trees)
{
  curline = taxamap[i,2:row_len]
  n_seqs = max(curline)
  curtree = trees[i]
  for (j in 2:n_seqs)
  {
    if (sum(curline == j) > 1)
    {
      dups = which(curline == j)
      curname = taxanames[dups[1]]
      newname = paste(taxanames[dups], collapse=",")
      stopifnot(length(grep(curname, curtree))>0)
      curtree = gsub(curname, newname, curtree)
    }
  }
  trees[i] = curtree

  perc_iter = perc_iter-1
  if (perc_iter <= 0)  
  {
    perc_count = perc_count+1
    perc_iter = perc_step
    cat(perc_count,"%\r")
  }
}
cat("... Done!                      \n")

write(trees, outputfile)

