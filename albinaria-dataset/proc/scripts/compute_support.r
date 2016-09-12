
fname = function(...)
{
  fn = paste(..., sep="")
  fn
}
rf_prefix = "locus.rfdists.mre."
nt_prefix = "perlocus.trees."
desc_file = "loci.desc"

####################################

locus_desc = scan(desc_file, list(locus_id=0, ntaxa=0, taxaprop=0, len=0, nvar=0, ninf=0, varprop=0, gapy=0, taxamap=""), comment.char="#", quiet=TRUE)

results = matrix(ncol=7, nrow=12046)
colnames(results) = c("id", "ntaxa", "ninf", "ntrees", "unique", "rf", "rsupport")

sequential_index = 1
for (klocus_id in 0:17)
{
  cat("parse file", klocus_id, "/ 17\n")
  rf_file    = fname(rf_prefix, klocus_id)
  ntrees_file = fname(nt_prefix, klocus_id)

  ntrees  = scan(ntrees_file, list(locus_id=0, total=0, unique=0, rf=0), quiet=TRUE)
  rfdists = scan(rf_file, list(tree_id=0, rf=0), quiet=TRUE)

  stopifnot(sum(ntrees$total) == length(rfdists$tree_id))

  n_locus = length(ntrees$locus_id)

  cur_id = 1
  for ( locus_id in 1:n_locus)
  {
    subset_start = cur_id
    subset_end   = cur_id + ntrees$total[locus_id] - 1

    r_locus_id = ntrees$locus_id[locus_id]
    r_ntaxa = locus_desc$ntaxa[r_locus_id]
    r_ninf = locus_desc$ninf[r_locus_id]
    r_ntrees = ntrees$total[locus_id]
    r_unique = ntrees$unique[locus_id]
    r_rf = mean(rfdists$rf[subset_start:subset_end])
    r_support = sum(rfdists$rf[subset_start:subset_end] == 0) / r_ntrees
    results[sequential_index,] = c(
        r_locus_id,
        r_ntaxa,
        r_ninf,
        r_ntrees, 
        r_unique, 
        r_rf,
        r_support) # % of trees with RF distance 0

    cur_id = subset_end + 1
    sequential_index = sequential_index + 1
  }
}

sig_results = results[results[,'ninf'] >= (results[,'ntaxa']-3),]

sum(sig_results[sig_results[,'unique'] == 1, 'rf'] == 0)

# avg number of informative sites
mean(sig_results[,'ninf'])
mean(sig_results[sig_results[,'unique']==1,'ninf'])
