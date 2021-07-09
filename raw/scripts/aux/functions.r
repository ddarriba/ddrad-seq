# Aux functions for R scripts

parse_filter = function(fname, data)
{
  values = scan(fname, 
                list(field="", value=0), 
                blank.lines.skip=TRUE,
                comment.char="#",
                quiet=TRUE)
  filter_cfg = list()
  filter_cfg$min_taxa = values$value[which(values$field == "min_taxa")]
  filter_cfg$min_var  = values$value[which(values$field == "min_var")]
  filter_cfg$min_inf  = values$value[which(values$field == "min_inf")]

  # build filter by taxa
  filter_cfg$taxa = data$eftaxa >= filter_cfg$min_taxa

  # build filter by number of variable sites
  if (filter_cfg$min_var < 0)
    filter_cfg$nvar = data$nvar   >= 1
  else
    filter_cfg$nvar = data$nvar   >= filter_cfg$min_var

  # build filter by number of informative sites
  if (filter_cfg$min_inf < 0)
    filter_cfg$ninf = data$ninf   >= 1
  else
    filter_cfg$ninf = data$ninf   >= filter_cfg$min_inf

  # build filter
  filter_cfg$all  = filter_cfg$taxa * filter_cfg$nvar * filter_cfg$ninf

  #print(sum(filter_cfg$taxa))
  #print(sum(filter_cfg$nvar))
  #print(sum(filter_cfg$ninf))
  #print(sum(filter_cfg$all))

  filter_cfg
}

parse_loci_desc = function(fname)
{
  loci_desc = scan(fname, 
                   list(id=0, ntax=0, tprop=0, len=0, 
                        nvar=0, ninf=0, vprop=0, gapy=0, 
                        tmap="", eftaxa=0, dups=""), 
                   comment.char="#", 
                   quiet=TRUE)

  loci_desc
}
