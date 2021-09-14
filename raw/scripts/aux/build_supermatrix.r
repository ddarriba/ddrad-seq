#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("scripts/aux/functions.r")

FORMAT = 0 # 0 for PHYLIP, 1 for FASTA

stopifnot(length(args) == 7)

#input files
base_name      = args[1]
data_dir       = args[2]
file_loci_desc = args[3]
file_taxa      = args[4]
dir_loci       = args[5]
file_filter    = args[6]
filter_name    = args[7]

#output data
output_base_dir = paste("../proc", data_dir, filter_name, sep='/')
file_parts = "supermatrix.parts"
if (FORMAT == 0) {
  file_msa = "supermatrix.phy"
} else {
  file_msa = "supermatrix.phy"
}
output_file_msa = paste(output_base_dir, file_msa, sep='/')
output_file_parts = paste(output_base_dir, file_parts, sep='/')

cat("Base name:      ", base_name, "\n")
cat("Data directory: ", data_dir, "\n")
cat("Loci desc file: ", file_loci_desc, "\n")
cat("Taxa file:      ", file_taxa, "\n")
cat("Loci directory: ", dir_loci, "\n")
cat("Filter file:    ", file_filter, "\n")
cat("Output:         ", output_file_msa, "\n")

tnames = scan(file_taxa, "", quiet=TRUE)

if (file.exists(output_file_msa) || file.exists(output_file_parts))
{
  unlink(output_file_msa)
  unlink(output_file_parts)
#  stop("Output files already exist: ", output_file_msa, ", ", output_file_parts, call.=FALSE)
}

#  id ntax  tprop   len  nvar   ninf  vprop   gapy   tmap eftaxa   dups
data = scan(file_loci_desc, list(id=0, ntax=0, tprop=0, len=0, nvar=0, ninf=0, vprop=0, gapy=0, tmap="", eftaxa=0, dups=""), comment.char="#", quiet=TRUE)

n_loci = length(data$id)
cat("Loci file contains", n_loci,"lines\n")

# apply filters
filter = parse_filter(file_filter, data)

cat("Original number of loci =", n_loci, "\n")
cat("Loci with less than", filter$min_taxa, "taxa =", n_loci - sum(filter$taxa), "\n")
cat("Loci with less than", filter$min_var, "variable sites =", n_loci - sum(filter$nvar), "\n")
cat("Loci with less than", filter$min_inf, "informative sites =", n_loci - sum(filter$ninf), "\n")
cat("Filtered number of loci =", sum(filter$all), "\n")

cat("\n\n")

n_alltaxa = length(strsplit(data$tmap[1], '')[[1]])

cat("Initializing supermatrix\n")
supermat_len = sum(data$len[filter$all == 1])
#supermat = matrix("?", ncol=supermat_len, nrow=n_alltaxa)
#rownames(supermat) = tnames

cat("Start processing\n")

if (FORMAT == 0) {
  write(x = paste(n_alltaxa, supermat_len), file = output_file_msa, append = FALSE)
}

# TODO AVOID ALLOCATING ALL MEMORY TOGETHER!!!

for (j in 1:n_alltaxa)
{
  cat("... progress T =",j, "/", n_alltaxa, "\n")
  curname = tnames[j]
  fastaname = paste(">",curname,sep='')
  k_id = 0
  msa_dir = paste(getwd(), dir_loci, 0, sep="/")
  locus_ini = 1
  curseq = matrix("?", ncol=supermat_len, nrow=1)
  for (id in 1:n_loci)
  {
    if (!(id%%1000))
    {
      k_id = k_id+1
      msa_dir = paste(getwd(), dir_loci, k_id, sep="/")
    }


    if (filter$all[id])
    {
      # process locus

      msa_filename = paste(base_name, "locus", id, sep=".")
      msa_file = paste(msa_dir, msa_filename, sep="/")
      stopifnot(file.exists(msa_file))

      locus_end = locus_ini + data$len[id] - 1

      if (j == 1)
      {
        # write partition
        write(x = paste("LOC",id," = ", locus_ini,"-",locus_end, sep=""), 
              file = output_file_parts, 
              append = TRUE)
      }

      msa_data = scan(msa_file, c(""), quiet=TRUE)
      msa_seqs = length(msa_data)/2
      stopifnot(msa_seqs == data$ntax[id])

      if (fastaname %in% msa_data[1:msa_seqs*2-1])
      {
        locus_index = which(msa_data[1:msa_seqs*2-1] == fastaname)*2
        seq = strsplit(msa_data[locus_index], '')[[1]]
        curseq[1, locus_ini:locus_end] = seq    
      }
      locus_ini = locus_end+1
    }
  }

  # Dump sequence
  if (FORMAT == 0) {
    write(x = paste(n_alltaxa, supermat_len), file = output_file_msa, append = FALSE)
    write(
       x     = paste(curname, paste(curseq[1,], collapse="")),
       file  = output_file_msa,
       append = TRUE)
  } else {
      write(
        x     = paste(">", curname, sep=""),
        file  = output_file_msa,
        append = TRUE)
      write(
        x     = paste(curseq[1,], collapse=""),
        file  = output_file_msa,
        append = TRUE)
  }
}

wrn=warnings()
if (!is.null(wrn))
  print(wrn)
