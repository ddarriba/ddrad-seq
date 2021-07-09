RAXML_BIN=raxmlHPC-SSE3
treesfile=$1
basedir=`echo $treesfile | rev | cut -d'/' -f2- | rev`
treesfname=`echo $treesfile | rev | cut -d'/' -f1 | rev`

logdir=$basedir/log
resdir=$basedir/res
mkdir -p $logdir $resdir

constype=STRICT
execname=$treesfname.$constype
$RAXML_BIN -m GTRGAMMA -J $constype -z $treesfile -n $execname
mv RAxML_info.$execname $logdir/$execname.log
mv RAxML_StrictConsensusTree.$execname $resdir/$execname

constype=MR
execname=$treesfname.$constype
$RAXML_BIN -m GTRGAMMA -J $constype -z $treesfile -n $execname
mv RAxML_info.$execname $logdir/$execname.log
mv RAxML_MajorityRuleConsensusTree.$execname $resdir/$execname

constype=MRE
execname=$treesfname.$constype
$RAXML_BIN -m GTRGAMMA -J $constype -z $treesfile -n $execname
mv RAxML_info.$execname $logdir/$execname.log
mv RAxML_MajorityRuleExtendedConsensusTree.$execname $resdir/$execname
