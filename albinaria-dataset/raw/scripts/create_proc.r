#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("scripts/functions.r")

#config
USE_SYMLINK     = TRUE  # whether to use a symbolic link or an actual copy
OVERRIDE_OUTPUT = TRUE

#input files
descfile = "loci.desc"
taxafile = "taxa.desc"
loci_dir = "loci"
filtercfgfile="scripts/filter.cfg"

#output data
output_base_dir = "../proc"
output_loci_dir = "loci"

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("Files prefix must be supplied", call.=FALSE)
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

output_dir=paste(output_base_dir, output_loci_dir, sep="/")
if (!dir.exists(output_base_dir)) dir.create(output_base_dir)
if (!dir.exists(output_dir)) dir.create(output_dir)

cat("Processing MSA files: output to",output_dir,"\n")
k_id=0
msa_dir=paste(getwd(), loci_dir, k_id, sep="/")
output_msa_dir=paste(output_dir, k_id, sep="/")
cat("MSA DIR =",msa_dir, "\n")
if (!dir.exists(output_msa_dir)) dir.create(output_msa_dir)

for (id in 1:n_loci)
{
  if (!(id%%1000))
  {
    cat("... progress K =",k_id+1, "/", n_loci%/%1000, "\n")
    k_id = k_id+1
    msa_dir = paste(getwd(), loci_dir, k_id, sep="/")
    output_msa_dir=paste(output_dir, k_id, sep="/")
    if (!dir.exists(output_msa_dir)) dir.create(output_msa_dir)
  }

  if (filter$all[id])
  {
    msa_filename = paste(prefix, "locus", id, sep=".")
    output_msa_file = paste(output_msa_dir, msa_filename, sep="/")

    if (data$ntax[id] != data$eftaxa[id])
    {
      msa_filename = paste(msa_filename,"reduced", sep=".")
    }

    msa_file = paste(msa_dir, msa_filename, sep="/")
    stopifnot(file.exists(msa_file))

    if (OVERRIDE_OUTPUT)
      unlink(output_msa_file)

    if (USE_SYMLINK)
      file.symlink(from=msa_file, to=output_msa_file)
    else
      file.copy(from=msa_file, to=output_msa_file)
  }
}

wrn=warnings()
if (!is.null(wrn))
  print(wrn)
