library(phangorn)
descfile    = "perlocus.trees"
taxamapfile = "taxamap"
taxafile    = "../raw/taxa"
dupsfile    = "taxa.dups"

outputfile = "loci.splitsupport"
outputtrees = "loci.toptrees"

do_output  = TRUE

if (do_output) {
  unlink(outputfile)
  unlink(outputtrees)
}

get_split = function(spl, indices, ntax, alltax, duplist)
{
  com_spl = numeric(alltax)
  com_spl[indices] = 1
  if (duplist[1] > 0)
    com_spl[duplist] = 1
  test_ind = 128
  for (i in 1:ntax)
  {
    k = (i-1)%/%8 + 1
#    if (!i%%8)
#     k2 = 0
#    else
#     k2 = 8 - i%%8
   
    exist_taxon = bitwAnd(spl[k], test_ind)
    if(exist_taxon)
    {
      com_spl[indices[i]] = 2
      if (indices[i] %in% duplist)
        com_spl[duplist[which(duplist == indices[i])+1]] = 2
    }
    test_ind = test_ind / 2
    if (test_ind < 1) test_ind = 128
  }
  com_spl
}

taxanames = scan(taxafile, c(""))
taxadups  = scan(dupsfile, c(""))
desc = scan(descfile, list(id=0, ntrees=0, ntopos=0, avg_rf=0), quiet=TRUE)
taxamap = read.table(file=taxamapfile)
n_loci = length(desc$id)
n_alltaxa = length(taxanames)

stopifnot(dim(taxamap)[1] == n_loci)
stopifnot(sum(desc$id != taxamap[,1]) == 0)

taxadups = taxadups[desc$id]

# exclude id column
taxamap = taxamap[,2:dim(taxamap)[2]]

# read tree
if (!file.exists(outputfile))
{
  for (i in 1:n_loci)
  {
    id = desc$id[i]
    k_id = id %/% 1000
    treesfile = paste("res/",k_id,"/Albinaria98inds_c91d6m4p3.locus.",id,".trees",sep="")
    trees = read.tree(treesfile)
    n_trees = length(trees)
    stopifnot(n_trees == desc$ntrees[i])
    n_topos = desc$ntopos[i]
    n_taxa = length(trees[[1]]$tip.label)

    if (n_topos == 1) {
      # dump tree
      write.tree(trees[[1]], file=outputtrees, append=TRUE)
    }

    # get all different splits + freqs
    splits = bitsplits(trees)
    n_splits = ncol(splits$matsplit)
    taxa_indices = sapply(splits$labels, function(x) which(taxanames == x))
    splits_table = matrix(nrow=n_splits, ncol=n_alltaxa)
    splits_support = splits$freq / n_trees

    duplist = as.numeric(strsplit(taxadups[i], split=",")[[1]])
    duptaxa = 0
    if (duplist[1] > 0)
      duptaxa = length(duplist)/2
    for (j in 1:n_splits)
    {
      splits_table[j,] = get_split(as.numeric(splits$matsplit[,j]), taxa_indices, n_taxa, n_alltaxa, duplist)
      write(file=outputfile, c(id, splits_support[j], paste(splits_table[j,], collapse="")), ncol=n_alltaxa+2, append=TRUE)
      stopifnot(sum(splits_table[j,]>0) == n_taxa + duptaxa)
    }
    print(i)
  }
}

# reload splits_table
splits = read.table("loci.splitsupport")
unlink("loci.splitsupport.user")
for (i in 1:(dim(splits)[1]))
{
  l=splits[i,3:100]
  l[l==0]="x"
  l[l==1]="0"
  l[l==2]="1"
  l=paste(l, collapse="")
  write(file="loci.splitsupport.user", c(splits[i,1], splits[i,2], l), ncol=3, append=TRUE)
}

hist(splits[,2], xlab="Support", ylab="Frequency", main="Split support")
