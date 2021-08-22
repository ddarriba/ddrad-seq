STATE_SCORES = c("A"=0,"C"=0,"G"=0,"T"=0,"U"=0,
                 "W"=2,"S"=2,"M"=2,"K"=2,"R"=2,"Y"=2,
                 "B"=4,"D"=4,"H"=4,"V"=4,"N"=16,"-"=16)
DEBUG = FALSE

strhash <- function(str)
{
  hv = as.integer(charToRaw(str))
  h = 0
  for (i in 1:length(hv))
    h = h + hv[i]*i
  h
}

weightsites <- function(mat, state_scores)
{
  for (i in 1:length(STATE_SCORES))
  {
    mat[mat==names(STATE_SCORES)[i]] = STATE_SCORES[i]
  }
  mat
}

# callback score functions
# 1. score sites according to ambiguities and missing data
snp.score.all <- function(locus)
{
  colSums(matrix(as.numeric(weightsites(locus)),ncol=ncol(locus)))
}
# 2. score sites according to missing data only
snp.score.nonmissing <- function(locus)
{
  colSums((locus == "N") + (locus == "-"))
}
# 3. simply get the first one
snp.score.first <- function(locus)
{
  scores = numeric(ncol(locus)) + 10
  scores[1] = 0

  scores
}

snp.score.nonmissing <- function(locus)
{
  scores = numeric(ncol(locus)) + 10
  scores[1] = 0

  scores
}

# select one snp per locus
snp.select <- function(loci_list, score_fun, informative_only=FALSE)
{
  snp_list = list()
  for (i in 1:loci_list$loci_count)
  {
     snp_list[[i]] = list()

     missing_data = NULL
     varsites = colnames(loci_list[[i]]$sequences)

     snp_list[[i]]$included = (loci_list[[i]]$var_count > 0) && (!informative_only || loci_list[[i]]$inf_count > 0)
     if (!snp_list[[i]]$included)
       next
       
     if (informative_only)     
       lmat = loci_list[[i]]$sequences[,varsites == "-I-"]
     else
       lmat = loci_list[[i]]$sequences[,varsites == "-I-" | varsites == "-V-"]

     stopifnot(!is.null(lmat))

     snp_list[[i]]$var_count = loci_list[[i]]$var_count
     snp_list[[i]]$inf_count = loci_list[[i]]$inf_count
     
     matdim = dim(lmat)
     stopifnot((loci_list[[i]]$var_count * loci_list[[i]]$taxa_count) == (matdim[1] * matdim[2]))
     
     if (is.null(dim(lmat)))
     {
       snp_list[[i]]$data = lmat
     }
     else
     {
       missing_data = score_fun(lmat)
       snp_list[[i]]$data = lmat[,which(missing_data == min(missing_data))[1]]
     }

     snp_list[[i]]$locus_id = loci_list[[i]]$locus_id     
     snp_list[[i]]$taxa_count = length(snp_list[[i]]$data)
     snp_list[[i]]$uniq_count = loci_list[[i]]$uniq_count
     stopifnot(snp_list[[i]]$taxa_count == loci_list[[i]]$taxa_count)
     stopifnot(snp_list[[i]]$taxa_count >= snp_list[[i]]$uniq_count)
  }
  snp_list$taxa = loci_list$taxa
  snp_list
}

# create snp matrix
snp.mat <- function(snp_list, snp_filter=NULL)
{
  snp_matrix = list()
  taxa_names = snp_list$taxa
  taxa_count = length(taxa_names)
  if (is.null(snp_filter))
    snp_filter = rep(TRUE, length(snp_list)-1)

  snp_count = sum(snp_filter)
  stopifnot(length(snp_filter) == (length(snp_list)-1))
  
  if (DEBUG)
    cat(snp_count, "SNPs included and", taxa_count, "taxa\n")
  
  snp_matrix$data = matrix(rep("-", taxa_count*snp_count), ncol=snp_count)
  rownames(snp_matrix$data) = taxa_names

  col_id = 1
  for (i in 1:length(snp_filter))
  {
    if (!is.null(snp_filter) && !snp_filter[i])
      next
    snp_matrix$data[names(snp_list[[i]]$data),col_id] = unlist(snp_list[[i]]$data)
    col_id = col_id+1
  }
  
  snp_matrix$parts = NULL
  snp_matrix
}

snp.filter <- function(loci_list, min_taxa, min_var, min_inf)
{
  loci_count   = loci_list$loci_count
  filter_array = rep(FALSE, loci_count)
  lmin_var = min_var
  lmin_inf = min_inf
  
  for (i in 1:loci_count)
  {
    if (min_var == -1) lmin_var = log2(loci_list[[i]]$taxa_count)
    if (min_inf == -1) lmin_inf = log2(loci_list[[i]]$taxa_count)
    filter_array[i] = loci_list[[i]]$taxa_count >= min_taxa && 
                      loci_list[[i]]$var_count  >= lmin_var && 
                      loci_list[[i]]$var_count  >= lmin_inf
  }
  
  filter_array
}

# create loci matrix
loci.mat <- function(loci_list, snp_filter=NULL)
{
  loci_matrix = list()
  
  taxa_names = loci_list$taxa
  taxa_count = length(taxa_names)
  if (is.null(snp_filter))
    snp_filter = rep(TRUE, loci_list$loci_count)

  snp_count = sum(snp_filter)
  stopifnot(length(snp_filter) == loci_list$loci_count)
  
  if (DEBUG)
    cat(snp_count, "SNPs included and", taxa_count, "taxa\n")

  matrix_length = 0
  for (i in 1:loci_list$loci_count)
    if (snp_filter[i])
      matrix_length = matrix_length + loci_list[[i]]$length
        
  loci_matrix$data = matrix(rep("-", taxa_count*matrix_length), ncol=matrix_length)
  loci_matrix$start = rep(0, sum(snp_filter))
  loci_matrix$end   = rep(0, sum(snp_filter))
  loci_matrix$id    = rep(0, sum(snp_filter))
  
  rownames(loci_matrix$data) = taxa_names

  col_id = 1
  next_part = 1
  for (i in 1:loci_list$loci_count)
  {
    if (!snp_filter[i])
      next
    len = loci_list[[i]]$length

    loci_matrix$data[rownames(loci_list[[i]]$sequences),col_id:(col_id+len-1)] = loci_list[[i]]$sequences
    loci_matrix$start[next_part] = col_id
    loci_matrix$end[next_part] = col_id + len - 1
    loci_matrix$id[next_part] = i
    
    next_part = next_part + 1
    col_id = col_id+len
  }
  loci_matrix
}

write.msa <- function(snp_matrix, output_file, format, snp_filter=NULL)
{
  if (is.null(snp_filter))
  {
    if (is.null(snp_matrix$start)
      snp_filter = rep(TRUE, dim(snp_list)[2])
    else
      snp_filter = rep(TRUE, length(snp_list$start))
  }
  
  if (tolower(format) == "phylip")
  {
    write(x = paste(dim(snp_matrix$data), collapse=" "), file = output_file, append=FALSE)
    for (i in 1:dim(snp_matrix$data)[1])
      write(
        x      = paste(rownames(snp_matrix$data)[i], paste(snp_matrix$data[i,], collapse="")),
        file   = output_file,
        append = TRUE)
      cat("MSA dumped to",output_file,"in PHYLIP format\n")
  }
  else if (tolower(format) == "fasta")
  {
    if (file.exists(output_file)) {
      file.remove(output_file)
    }
    for (i in 1:dim(snp_matrix$data)[1])
    {
      write(
        x      = paste(">", rownames(snp_matrix$data)[i], sep=""),
        file   = output_file,
        append = TRUE)
      write(
        x      = paste(snp_matrix$data[i,], collapse=""),
        file   = output_file,
        append = TRUE)
    }
    cat("MSA dumped to",output_file,"in FASTA format\n")
  }
  else
  {
    cat("Unknown format",format,". It should be FASTA or PHYLIP\n")
    return
  }
  
  if (!is.null(snp_matrix$start))
  {
    parts_file = paste(output_file,"parts",sep=".")
    # dump partitions file
    if (file.exists(parts_file)) {
      file.remove(parts_file)
    }
    for (i in 1:length(snp_matrix$start))
      write(
        x      = paste("PART", snp_matrix$id[i],"=",snp_matrix$start[i],"-",snp_matrix$end[i], sep=""),
        file   = parts_file,
        append = TRUE)
      cat("Partitions dumped to", parts_file,"in RAxML format\n")
  }
}


read.loci <- function(filepath, collapse_identical_sequences=TRUE)
{
  loci_file = file(filepath, "r")
  on.exit(close(loci_file))
  
  # global
  loci_data = NULL
  loci_all_count = 0
  loci_inc_count = 0
  gtaxa = c()
  
  # local
  sequences = c()
  lhash = c()
  ltaxa = c()
  lseq_count = 0
  luniq_count = 0
  include_locus = TRUE  
  locus_length = 0
  seq_start = 0
  
  while(TRUE) {
    line = readLines(loci_file, n=1)
    if ( length(line) == 0 )
      break
      
    if (grepl('^>', line))
    {
      # sequence
      lseq_count = lseq_count + 1
            
      pieces=strsplit(line, "[ |\t]+")
      ltaxa[lseq_count] = substring(pieces[[1]][1], first=2)

      sequences[lseq_count] = pieces[[1]][2]
      lh = strhash(sequences[lseq_count])
      
      if (lseq_count == 1)
      {
         # first line
         locus_length = nchar(sequences[lseq_count])
         seq_start = gregexpr(pattern=sequences[lseq_count],line)[[1]][1]
         line_end = nchar(line)
         lhash[1] = lh
         luniq_count = 1
      }
      else
      {
         dup = FALSE
         for (i in 1:length(lhash))
         {
           if (lhash[i] == lh)
           {
             if (sequences[lseq_count] == sequences[i])
             {
               dup = TRUE
               break
             }
           }
         }
         lhash[lseq_count] = lh
         
         if (!dup)
           luniq_count = luniq_count + 1
      }
    }
    else  if (grepl('^//', line))
    {
      # meta

      # .. correct 'meta' start
      match_end = gregexpr(pattern="\\|",line)[[1]][1]
      seq_start = seq_start + (line_end-match_end-1)

      metainfo = substring(line, first=seq_start, last=match_end-1)
      
      locus_id = as.integer(strsplit(line, "\\|")[[1]][2])
      
      loci_all_count = loci_all_count + 1
      if (loci_all_count %% 1000 == 0) print(loci_all_count)

      inf_sites=gregexpr(pattern='\\*', metainfo)[[1]]
      if (inf_sites[1] > -1)
        inf_count = length(inf_sites)
      else
        inf_count = 0
      var_sites=gregexpr(pattern='\\-', metainfo)[[1]]
      if (var_sites[1] > -1)
        var_count = length(var_sites)
      else
        var_count = 0
      var_count = var_count + inf_count

      # .. compose entry
      loci_data[[loci_all_count]] = list()
      loci_data[[loci_all_count]]$locus_id = locus_id
      loci_data[[loci_all_count]]$length = locus_length
      loci_data[[loci_all_count]]$taxa_count = length(ltaxa)
      loci_data[[loci_all_count]]$uniq_count = luniq_count
      loci_data[[loci_all_count]]$var_sites = var_sites
      loci_data[[loci_all_count]]$var_count = var_count
      loci_data[[loci_all_count]]$inf_sites = inf_sites
      loci_data[[loci_all_count]]$inf_count = inf_count
      #loci_data[[loci_all_count]]$sequences = sequences
      loci_data[[loci_all_count]]$sequences = matrix(strsplit(paste(sequences, collapse=""), "")[[1]],ncol=locus_length,byrow=TRUE)
      rownames(loci_data[[loci_all_count]]$sequences) = ltaxa

      cn = rep("", locus_length)
      if (inf_count > 0)
        cn[inf_sites] = "-I-"
      if (var_count > inf_count)
        cn[var_sites] = "-V-"

      colnames(loci_data[[loci_all_count]]$sequences) = cn
      
      gtaxa = unique(c(gtaxa, ltaxa))
      
      # .. reset
      sequences = c()
      ltaxa = c()
      lseq_count = 0
      luniq_count = 0
      lhash = c()
      include_locus = TRUE
    }
  }

  cat("Read", loci_all_count, "loci\n")
  if (DEBUG)
    cat("Close File\n")
  
  # return
  loci_data$taxa = gtaxa
  loci_data$loci_count = loci_all_count
  loci_data
}

