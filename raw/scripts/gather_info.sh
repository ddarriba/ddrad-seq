#!/bin/bash
# This script gathers per-locus info from .loci file
# 
# %.loci.head file contains the headers of %.loci but taxa names are
# converted to the indices in 'taxa' file.
#
# Run: gather_info.sh [DATA_DIR]
# 

source scripts/aux/common.sh

RAXML_EXECNAME=readseq
RAXML_file_loci_desc=RAxML_info.$RAXML_EXECNAME

# input
[[ -f ${file_loci} ]] || { echo "Loci file missing"; exit; }
[[ -f ${file_loci_head} ]] || { echo "Loci head file missing"; exit; }
[[ -f ${file_taxa} ]] || { echo "Taxa file missing"; exit; }

# output
[[ -z ${dir_loci} ]] && { echo "Loci directory undefined"; exit; }
[[ -z ${file_loci_desc} ]] && { echo "Loci desc file undefined"; exit; }

mkdir -p $dir_loci

max_taxa=`cat ${file_taxa} | wc -l`

header_lines_count=13 # count of lines printed above

loci_sections=`fgrep -n "//" $file_loci | cut -d':' -f 1`
cur_line=1
locus_id=1
n_loci=`echo $loci_sections | wc -w`

if [ -f ${file_loci_desc} ]; then
  # check if file is complete
  loci_desc_lines=`cat ${file_loci_desc} | wc -l`
  [[ $((loci_desc_lines - header_lines_count)) == ${n_loci} ]] && \
    { echo "Loci desc file already exist"; exit; } || \
    { echo "Loci desc file contains wront number of lines"; }
fi

echo "# id:     locus index" > $file_loci_desc
echo "# ntax:   number of taxa" >> $file_loci_desc
echo "# tprop:  proportion of taxa (ntax/98)" >> $file_loci_desc
echo "# len:    length (number of sites)" >> $file_loci_desc
echo "# nvar:   number of variable sites" >> $file_loci_desc
echo "# ninf:   number of informative sites" >> $file_loci_desc
echo "# nbsnp:  number of biallelic snps" >> $file_loci_desc
echo "# vprop:  proportion of var sites (nvar + ninf)/len" >> $file_loci_desc
echo "# gapy:   gapyness (proportion of gaps: n_gaps/(ntax*len))" >> $file_loci_desc
echo "# tmap:   existing taxa in the locus" >> $file_loci_desc
echo "# eftax:  effective number of taxa" >> $file_loci_desc
echo "# dups:   list of duplicated sequences" >> $file_loci_desc
printf "#%4s %4s %6s %5s %5s %6s %6s %6s %6s %6s %6s\n" id ntax tprop len nvar ninf vprop gapy tmap eftaxa dups >> $file_loci_desc

for locus_end in $loci_sections; do
  locus_file=${dir_loci}/$((locus_id/1000))/${base_name}.locus.${locus_id}

  [[ -f ${locus_file} ]] || { echo "Locus file ${locus_file} does not exist"; exit; }

  # process locus_file with raxml
  $RAXML_BIN -f c -s ${locus_file} -m GTRGAMMA -n $RAXML_EXECNAME > /dev/null

  [[ -f ${RAXML_file_loci_desc} ]] || { echo "RAxML output ${RAXML_file_loci_desc} does not exist"; exit; }

  identical_seqs=`grep "WARNING:" $RAXML_file_loci_desc | grep "exactly identical" | cut -d' ' -f4,6`
  ndups=`echo $identical_seqs | wc -w`
  dups="0,0"
  for tname in $identical_seqs; do
    tid=`grep -n "^${tname}$" $file_taxa | cut -d':' -f 1`
    dups="$dups,$tid"
  done

  rm $RAXML_file_loci_desc

  n_taxa=$((locus_end-cur_line))
  eftaxa=$((n_taxa - ndups/2))

  sed -n "${cur_line},$((locus_end))p" $file_loci | rev > $file_tmp

  n_sites=`sed "${cur_line}q;d" ${file_loci} | xargs | cut -d' ' -f2 | wc -c`
  n_sites=$((n_sites - 1))
  varsites_line=`sed "${locus_end}q;d" ${file_loci}`
  n_var=$((`echo "${varsites_line//[^'-']}" | wc -c`-1))
  n_inf=$((`echo "${varsites_line//[^'*']}" | wc -c`-1))
  var_prop=`echo "scale=4; ($n_var + $n_inf)/$n_sites" | bc -l`
  taxa_prop=`echo "scale=4; $n_taxa/$max_taxa" | bc -l`

  taxa_inc=`sed -n ${cur_line},$((locus_end-1))p ${file_loci_head} | sort -g`
  
  #TODO: check that words in taxa_inc equals number of taxa
  cnt=1
  tmap=
  for i in $taxa_inc; do
    for j in `seq $((cnt+1)) $i`; do
      tmap=${tmap}0
    done
    cnt=$((i+1))
    tmap=${tmap}1
  done
  for j in `seq $cnt $max_taxa`; do
    tmap=${tmap}0
  done

  n_gaps=`fgrep -o '-' $locus_file | wc -l`
  gapyness=`echo "scale=4; $n_gaps/($n_taxa * $n_sites)" | bc -l`

  printf "%5s %4s %6s %5s %5s %6s %6s %6s   %s %6s %s\n" $locus_id $n_taxa $taxa_prop $n_sites $n_var $n_inf $var_prop $gapyness $tmap $eftaxa $dups >> $file_loci_desc
  echo ${locus_id}/${n_loci}

  locus_id=$((locus_id+1))
  cur_line=$((locus_end+1))
done
