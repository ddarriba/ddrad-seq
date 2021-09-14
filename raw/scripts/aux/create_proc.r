#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("scripts/aux/functions.r")

#config
USE_SYMLINK     = TRUE  # whether to use a symbolic link or an actual copy
OVERRIDE_OUTPUT = TRUE

stopifnot(length(args) == 6)

# input files
base_name       = args[1]
data_dir        = args[2]
file_loci_desc  = args[3]
dir_loci        = args[4]
file_filter_cfg = args[5]
filter_name     = args[6]

#output data
output_base_dir = paste("../proc", data_dir, filter_name, sep='/')
output_loci_dir = paste(output_base_dir, "loci", sep='/')

#  id ntax  tprop   len  nvar   ninf  vprop   gapy   tmap eftaxa   dups
data = scan(file_loci_desc, list(id=0, ntax=0, tprop=0, len=0, nvar=0, ninf=0, vprop=0, gapy=0, tmap="", eftaxa=0, dups=""), comment.char="#", quiet=TRUE)

n_loci = length(data$id)
cat("Loci file contains", n_loci,"lines\n")

# apply filters
filter = parse_filter(file_filter_cfg, data)

cat("Original number of loci =", n_loci, "\n")
cat("Loci with less than", filter$min_taxa, "taxa =", n_loci - sum(filter$taxa), "\n")
cat("Loci with less than", filter$min_var, "variable sites =", n_loci - sum(filter$nvar), "\n")
cat("Loci with less than", filter$min_inf, "informative sites =", n_loci - sum(filter$ninf), "\n")
cat("Filtered number of loci =", sum(filter$all), "\n")

cat("Output directory =", output_base_dir, "\n")
cat("Output loci directory =", output_loci_dir, "\n")

cat("\n\n")

if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive=TRUE)
if (!dir.exists(output_loci_dir)) dir.create(output_loci_dir, recursive=TRUE)

cat("Processing MSA files: output to",output_loci_dir,"\n")
k_id=0
msa_dir=paste(getwd(), dir_loci, k_id, sep="/")
output_msa_dir=paste(output_loci_dir, k_id, sep="/")

cat("MSA DIR =",msa_dir, "\n")

if (!dir.exists(output_msa_dir)) dir.create(output_msa_dir)

for (id in 1:n_loci)
{
  if (!(id%%1000))
  {
    k_id = k_id+1
    msa_dir = paste(getwd(), dir_loci, k_id, sep="/")
    output_msa_dir=paste(output_loci_dir, k_id, sep="/")
    if (!dir.exists(output_msa_dir)) dir.create(output_msa_dir)
    cat("... progress K =",k_id+1, "/", n_loci%/%1000, " : ", output_msa_dir, "\n")
  }

  if (filter$all[id])
  {
    msa_filename = paste(base_name, "locus", id, sep=".")
    output_msa_file = paste(output_msa_dir, msa_filename, sep="/")

    cat(msa_filename, data$id[id], data$ntax[id], data$eftax[id], "\n")
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
