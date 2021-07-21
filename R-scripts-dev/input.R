STATE_SCORES = c("A"=0,"C"=0,"G"=0,"T"=0,"U"=0,
                 "W"=2,"S"=2,"M"=2,"K"=2,"R"=2,"Y"=2,
                 "B"=4,"D"=4,"H"=4,"V"=4,"N"=16,"-"=16)
DEBUG = FALSE

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
snp.select <- function(loci_list, score_fun)
{
  snp_list = list()
  for (i in 1:length(loci_list))
  {
     missing_data = NULL
     lmat = loci_list[[i]]
     if (is.null(lmat)) next

          
     if (is.null(dim(lmat)))
     {
       snp_list[[i]] = lmat
     }
     else
     {
       missing_data = score_fun(lmat)
       snp_list[[i]] = lmat[,which(missing_data == min(missing_data))[1]]
     }
  }
  snp_list
}

# create snp matrix
snp.mat <- function(snp_list, taxa_filepath, min_taxa=3)
{
  taxa_names = scan(taxa_filepath, "", quiet=TRUE)
  taxa_count = length(taxa_names)
  snp_count = sum(unlist(lapply(snp_list, function(x) length(x)>=min_taxa)))
  
  if (DEBUG)
    cat(snp_count, "SNPs found and", taxa_count, "taxa\n")
  
  snp_matrix = matrix(rep("-", taxa_count*snp_count), ncol=snp_count)
  rownames(snp_matrix) = taxa_names

  col_id = 1
  for (i in 1:length(snp_list))
  {
    if ((!is.null(snp_list[[i]])) && (length(snp_list[[i]]) >= min_taxa))
    {
      snp_matrix[names(snp_list[[i]]),col_id] = unlist(snp_list[[i]])
#      colnames(snp_matrix)[col_id] = i
      col_id = col_id+1
    }
  }
  snp_matrix
}

snp.dumpmatrix <- function(snp_matrix, output_file, format)
{
  if (tolower(format) == "phylip")
  {
    write(x = paste(dim(snp_matrix), collapse=" "), file = output_file, append=FALSE)
    for (i in 1:dim(snp_matrix)[1])
      write(
        x      = paste(rownames(snp_matrix)[i], paste(snp_matrix[i,], collapse="")),
        file   = output_file,
        append = TRUE)
      cat("MSA dumped to",output_file,"in PHYLIP format\n")
  }
  else if (tolower(format) == "fasta")
  {
    if (file.exists(output_file)) {
      file.remove(output_file)
    }
    for (i in 1:dim(snp_matrix)[1])
    {
      write(
        x      = paste(">", rownames(snp_matrix)[i], sep=""),
        file   = output_file,
        append = TRUE)
      write(
        x      = paste(snp_matrix[i,], collapse=""),
        file   = output_file,
        append = TRUE)
    }
    cat("MSA dumped to",output_file,"in FASTA format\n")
  }
  else
  {
    cat("Unknown format",format,". It should be FASTA or PHYLIP\n")
  }
}

read.loci <- function(filepath)
{
  loci_file = file(filepath, "r")
  taxa_count = 0
  seq_start = 0
  

  taxa = c()
  sequences = c()
  all_loci = NULL
  include_locus = TRUE
      
  while(TRUE) {
    line = readLines(loci_file, n=1)
    if ( length(line) == 0 )
      break
    
    if (grepl('^>', line))
    {
      pieces=strsplit(line, "[ |\t]+")
      taxon_name = substring(pieces[[1]][1], first=2)
      sequence = pieces[[1]][2]
      
      taxa = c(taxa, taxon_name)
      sequences = c(sequences,strsplit(sequence, "")[[1]])

      if (taxa_count == 0)
      {
         # first line
         loci_len = nchar(sequence)
         seq_start = gregexpr(pattern=sequence,line)[[1]][1]
      }
      taxa_count = taxa_count+1
    }
    else if (grepl('^//', line))
    {
      metainfo = substring(line, first=seq_start, last=seq_start+loci_len-1)
      loci_id = as.integer(substring(line, first=seq_start+loci_len+1))
      
      info_sites=gregexpr(pattern='\\*', metainfo)[[1]]

      if (info_sites[1] < 0)
      {
        # locus is not suitable
        info_sites_count = 0
        include_locus = FALSE
      }
      
      if (include_locus)
      {
        # locus is suitable
	info_sites_count = length(info_sites)

        loci = matrix(sequences, ncol=loci_len, byrow=TRUE)	
        rownames(loci) = taxa

        if (DEBUG)
        {
          print(loci[,info_sites])
          cat(loci_id, taxa_count, loci_len, nchar(metainfo), "INFO: ", info_sites_count,"\n")
          cat("include locus",loci_id,taxa_count,"\n")
        }
      }
      else
      {
        loci = NULL
        if (DEBUG)
          cat("skip locus",loci_id,"\n")
      }
      
      all_loci[[loci_id]] = loci[,info_sites]
 
      # reset
      taxa_count=0
      taxa = c()
      sequences = list()
      include_locus = TRUE
    }
    else
    {
      print("ERROR")
      break
    }
  }
  if (DEBUG)
    cat("Close File\n")
  close(loci_file)
  
  all_loci
}
