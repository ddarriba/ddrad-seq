#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("scripts/functions.r")

FORMAT = 0 # 0 for PHYLIP, 1 for FASTA
OUTPUT_MSA = "supermatrix"
OUTPUT_PARTS = "supermatrix.parts"

if (FORMAT == 0) {
  OUTPUT_MSA = paste(OUTPUT_MSA, "phy", sep=".")
} else {
  OUTPUT_MSA = paste(OUTPUT_MSA, "fas", sep=".")
}

#input files
descfile = "loci.desc"
taxafile = "taxa"
loci_dir = "loci"
filtercfgfile="scripts/filter.cfg"

tnames = scan(taxafile, "", quiet=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("Files prefix must be supplied", call.=FALSE)
}

if (file.exists(OUTPUT_MSA) || file.exists(OUTPUT_PARTS))
{
  stop("Output files already exist: ", OUTPUT_MSA, ", ", OUTPUT_PARTS, call.=FALSE)
}

prefix=args[1]

#  id ntax  tprop   len  nvar   ninf  vprop   gapy   tmap eftaxa   dups
data = scan(descfile, list(id=0, ntax=0, tprop=0, len=0, nvar=0, ninf=0, vprop=0, gapy=0, tmap="", eftaxa=0, dups=""), comment.char="#", quiet=TRUE)

n_loci = length(data$id)
cat("Loci file contains", n_loci,"lines\n")

# apply filters
filter = parse_filter(filtercfgfile, data)

cat("Original number of loci =", n_loci, "\n")
cat("Loci with less than", filter$min_taxa, "taxa =", n_loci - sum(filter$taxa), "\n")
cat("Loci with less than", filter$min_var, "variable sites =", n_loci - sum(filter$nvar), "\n")
cat("Loci with less than", filter$min_inf, "informative sites =", n_loci - sum(filter$ninf), "\n")
cat("Filtered number of loci =", sum(filter$all), "\n")

cat("\n\n")

n_alltaxa = length(strsplit(data$tmap[1], '')[[1]])

cat("Initializing supermatrix\n")
supermat_len = sum(data$len[filter$all == 1])
supermat = matrix("?", ncol=supermat_len, nrow=n_alltaxa)
rownames(supermat) = tnames

cat("Start processing\n")
k_id = 0
msa_dir = paste(getwd(), loci_dir, 0, sep="/")
locus_ini = 1
for (id in 1:n_loci)
{
  if (!(id%%1000))
  {
    cat("... progress K =",k_id+1, "/", n_loci%/%1000, "\n")
    k_id = k_id+1
    msa_dir = paste(getwd(), loci_dir, k_id, sep="/")
  }


  if (filter$all[id])
  {
    msa_filename = paste(prefix, "locus", id, sep=".")
    msa_file = paste(msa_dir, msa_filename, sep="/")
    stopifnot(file.exists(msa_file))

    locus_end = locus_ini + data$len[id] - 1

    write(x = paste("LOC",id," = ", locus_ini,"-",locus_end, sep=""), 
          file = OUTPUT_PARTS, 
          append = TRUE)

    msa_data = scan(msa_file, c(""), quiet=TRUE)
    msa_seqs = length(msa_data)/2
    stopifnot(msa_seqs == data$ntax[id])

    for (j in 1:msa_seqs)
    {
      tname = substring(msa_data[2*j - 1], 2)
      seq = strsplit(msa_data[2*j], '')[[1]]
      tindex = which(tnames == tname)
      supermat[tindex, locus_ini:locus_end] = seq  
    }
    locus_ini = locus_end+1
  }
}

cat("Dump MSA to",OUTPUT_MSA,"\n")
if (FORMAT == 0) {
  write(x = paste(n_alltaxa, supermat_len), file = OUTPUT_MSA, append = FALSE)
  for (taxon in 1:n_alltaxa)
  {
    cat("... sequence",taxon,"/",n_alltaxa,"\n")
    write(
      x     = paste(tnames[taxon], paste(supermat[taxon,], collapse="")),
      file  = OUTPUT_MSA,
      append = TRUE)
  }
} else {
  for (taxon in 1:n_alltaxa)
  {
    cat("... sequence",taxon,"/",n_alltaxa,"\n")
    write(
      x     = paste(">", tnames[taxon], sep=""),
      file  = OUTPUT_MSA,
      append = TRUE)
    write(
      x     = paste(supermat[taxon,], collapse=""),
      file  = OUTPUT_MSA,
      append = TRUE)
  }
}

wrn=warnings()
if (!is.null(wrn))
  print(wrn)
