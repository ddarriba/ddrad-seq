#!/bin/bash
# gather per-locus info from .loci file

base_name=Albinaria98inds_c91d6m4p3
loci_file=${base_name}.loci
perlocus_dir=perlocus

outfile=loci.desc
taxa_file=taxa
max_taxa=98

echo "# id:     locus index" > $outfile
echo "# ntax:   number of taxa" >> $outfile
echo "# tprop:  proportion of taxa (ntax/98)" >> $outfile
echo "# len:    length (number of sites)" >> $outfile
echo "# nvar:   number of variable sites" >> $outfile
echo "# ninf:   number of informative sites" >> $outfile
echo "# vprop:  proportion of var sites (nvar + ninf)/len" >> $outfile
echo "# gapy:   gapyness (proportion of gaps: n_gaps/(ntax*len))" >> $outfile
printf "#%4s %4s %6s %5s %5s %6s %6s %6s\n" id ntax tprop len nvar ninf vprop gapy >> $outfile

loci_sections=`fgrep -n "//" $loci_file | cut -d':' -f 1`
cur_line=1
locus_id=1
n_loci=`echo $loci_sections | wc -w`
for locus_end in $loci_sections; do

  n_taxa=$((locus_end-cur_line))
  n_sites=`sed "${cur_line}q;d" ${loci_file} | xargs | cut -d' ' -f2 | wc -c`
  varsites_line=`sed "${locus_end}q;d" ${loci_file}`
  n_var=$((`echo "${varsites_line//[^'-']}" | wc -c`-1))
  n_inf=$((`echo "${varsites_line//[^'*']}" | wc -c`-1))
  var_prop=`echo "scale=4; ($n_var + $n_inf)/$n_sites" | bc -l`
  taxa_prop=`echo "scale=4; $n_taxa/$max_taxa" | bc -l`

  locus_file=${perlocus_dir}/$((locus_id/1000))/${base_name}.locus.${locus_id}
  n_gaps=`fgrep -o '-' $locus_file | wc -l`
  gapyness=`echo "scale=4; $n_gaps/($n_taxa * $n_sites)" | bc -l`

  printf "%5s %4s %6s %5s %5s %6s %6s %6s\n" $locus_id $n_taxa $taxa_prop $n_sites $n_var $n_inf $var_prop $gapyness >> $outfile
  echo ${locus_id}/${n_loci}

  locus_id=$((locus_id+1))
  cur_line=$((locus_end+1))
done

exit
# find taxa representation
i=0
rm taxa.desc
while read taxon; do
  repr=`fgrep $taxon $loci_file | wc -l`
  rperc=`echo "scale=4; 100*$repr/$n_loci" | bc -l`
  printf "%12s %8s %8s\n" $taxon $repr $rperc >> taxa.desc
  echo $i
  i=$((i+1))
done < taxa
