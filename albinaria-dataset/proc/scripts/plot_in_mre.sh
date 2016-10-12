#!/bin/sh

# Plot all per-locus trees into the MRE tree
#
# Call:
#   $ plot_in_mre.sh [MRE_TREE]
#
# Requires:
#   * perlocus trees: res/x/y.trees
#
# Output:
#   * locus.rfdists.mre.%

RAXML_BIN=raxmlHPC-SSE3
mre_tree=$1
output_prefix="locus.rfdists.mre"

if [ -z $mre_tree ]; then
  echo Call plot_in_mre.sh [MRE_TREE]
  exit
fi

if [ ! -f $mre_tree ]; then
  echo MRE tree does not exist 
  exit
fi

kdirs=`ls res | wc -l`

for k_id in `ls res`; do
  start=$((k_id * 1000))
  end=$((k_id * 1000 + 999))
  rm -f tmpfile
  for i in `seq $start $end`; do
    f=res/$k_id/*.$i.trees
    [ `echo $f | wc -w` != 1 ] && echo ERROR && exit
    [ -f $f ] && cat $f >> tmpfile
  done

  execname=mre.rf.$k_id
  $RAXML_BIN -m GTRGAMMA -f R -t $mre_tree -z tmpfile -n $execname > /dev/null
  mv RAxML_RF-Distances.$execname ${output_prefix}.${k_id}

  echo $k_id
done

rm -f tmpfile
