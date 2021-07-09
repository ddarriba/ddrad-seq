#!/bin/sh

PREFIX=$1
kindex=$2
RES_DIR=$3
LOG_DIR=$4
OUT_DIR=$5

[ -z $OUT_DIR ] && exit

RAXML_BIN=raxmlHPC-SSE3
resdir=$RES_DIR/$kindex
logdir=$LOG_DIR/consensus/$kindex
outfile=$OUT_DIR/consensus.$kindex.trees
tmpfile=/tmp/tmp.lociconsensus.$kindex

CONSENSUS_TYPE=MRE
RAXML_RESULT_PREFIX=RAxML_MajorityRuleExtendedConsensusTree

mkdir -p $logdir
rm -f $outfile

for tfile in $resdir/$PREFIX.locus.*.trees; do
  i=`echo $tfile | rev | cut -d'.' -f2 | rev`
  echo "build consensus for $tfile : $i"
  execname=locusmre.$i
  $RAXML_BIN -m GTRGAMMA -J $CONSENSUS_TYPE -z $tfile -n $execname 2>&1
  tree=`cat $RAXML_RESULT_PREFIX.$execname`
  echo $i $tree >> $outfile
  mv RAxML_info.$execname $logdir
  mv $RAXML_RESULT_PREFIX.$execname $resdir
  echo "done $i"
done

mv $outfile $tmpfile
sort -n $tmpfile > $outfile
rm $tmpfile
