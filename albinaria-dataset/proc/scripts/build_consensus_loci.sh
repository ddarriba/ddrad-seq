kindex=$1

RAXML_BIN=raxmlHPC-SSE3
resdir=res/$kindex
logdir=log/consensus/$kindex
outfile=consensus.$kindex.trees

mkdir -p $logdir
rm -f $outfile

kstart=$((kindex * 1000))
kend=$(((kindex+1) * 1000 - 1))
for i in `seq $kstart $kend`; do
  tfile=`ls $resdir/*.locus.$i.trees 2> /dev/null`
  if [ ! -z $tfile -a -f $tfile ]; then
    execname=strict.$i
    $RAXML_BIN -m GTRGAMMA -J STRICT -z $tfile -n $execname 2>&1 > /dev/null
    tree=`cat RAxML_StrictConsensusTree.$execname`
    echo $i $tree >> $outfile
    mv RAxML_info.$execname $logdir
    mv RAxML_StrictConsensusTree.$execname $logdir
    echo $i
  fi
done
