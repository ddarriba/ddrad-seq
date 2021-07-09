#!/bin/bash

# Usage:
#
# ntrees.sh SUBDIR K_DIR
#
#    SUBDIR = data set subdirectory in data/
#    K_DIR = 0..(n_kdirs-1)
#

source scripts/aux/common.sh

raxml_bin="/usr/local/bin/raxmlHPC-SSE3"
outfile_base=perlocus.trees

n_kdirs=`ls res | wc -l`
kdirs=`seq 0 $((n_kdirs-1))`

if [[ ! -f "$raxml_bin" ]]; then
  echo "Error: Binary $raxml_bin does not exist."
  exit
fi

rm -f ${outfile_base}

for kvalue in $kdirs; do
  outfile=${outfile_base}.${kvalue}
  echo $kvalue : Dump into ${outfile}
  rm -f ${outfile}
  start_index=$((1000 * $kvalue))
  end_index=$((1000 * ($kvalue + 1)))

	resFolder=${dir_results}/$kvalue
  datFN=`echo $datFolder | rev | cut -d'/' -f1 | rev`
	for ((i=$start_index; i<$end_index; i++)); do #treesFile in $resFolder/*.trees; do
    treesFile="${resFolder}/${base_name}.locus.${i}.*"
    treesFile=`ls $treesFile 2> /dev/null`
    if [ ! -z $treesFile ]; then
      datIndex=$i #((kvalue*1000 + i))
      execname=ntrees.$datIndex
	  	$RAXML_BIN -f r -z $treesFile 0 -n $execname -m GTRGAMMA > /dev/null
      rax_outfile="RAxML_info.$execname"
      rffile="RAxML_RF-Distances.$execname"
      n_trees=`fgrep "Number of unique trees" $rax_outfile | rev | cut -d' ' -f 1 | rev`
      c_trees=`cat $treesFile | wc -l`
      values=`cat $rffile | cut -d' ' -f 3`
      pairs=`echo $values | wc -w`
      values=`echo $values | sed "s/ /+/g"`
      mean=`echo "scale=2; ($values)/$pairs" | bc -l`
      echo $datIndex $c_trees $n_trees $mean >> ${outfile}
      rm RAxML_*.$datIndex
    fi
	done
  cat ${outfile} >> ${outfile_base}
done
