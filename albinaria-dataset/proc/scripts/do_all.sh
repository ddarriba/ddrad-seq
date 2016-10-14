#!/bin/bash

RAW_DIR=../raw
LOCI_DIR=loci
RES_DIR=res
LOG_DIR=log
SM_DIR=supermatrices
CONSENSUS_DIR=loci_consensus
FILTERS_DIR=filters

TAXA_FILE=taxa
LOCIDESC_FILE=loci.desc

# used scripts
SCRIPT_EVAL_LOCI=scripts/eval_loci.sh
SCRIPT_BUILD_SM=scripts/build_supermatrix.r
SCRIPT_LOCI_CONSENSUS=scripts/build_consensus_loci.sh

function exit_error {
  echo ... error $1
  exit
}

function test_dir {
  # expect $1: directory name
  # expect $2: verbose (true | false)
  dirname=$1
  verbose=$2
  if [ ! -d $dirname ]; then
    $verbose && echo check directory \"$1\": NOPE
    return 1
  else
    $verbose && echo check directory \"$1\": OK!
    return 0
  fi
}

function test_file {
  # expect $1: file name
  # expect $2: verbose (true | false)
  filename=$1
  verbose=$2
  if [ ! -f $filename ]; then
    $verbose && echo check file \"$1\": NOPE
    return 1
  else
    $verbose && echo check file \"$1\": OK!
    return 0
  fi
}

function ask_confirm {
  # excpect $*: prompt
  prompt=$1
  read -p "$*" confirm
  [ -z $confirm ] && return 0
  [ "$confirm" == "Y" -o "$confirm" == "y" ] && return 0
  return 1
}

# step 0: check structure and get parameters

test_dir $LOCI_DIR true || exit_error
n_kdirs=`ls $LOCI_DIR | wc -l`
[ "$n_kdirs" == 0 ] && exit_error ": loci dir is empty"
for i in `seq 0 $((n_kdirs-1))`; do
  test_dir $LOCI_DIR/$i false || exit_error ": missing $LOCI_DIR/$i"
done
prefix=`ls $LOCI_DIR/0 | head -n 1 | cut -d'.' -f1`

test_file $TAXA_FILE true || \
  ( test_file $RAW_DIR/$TAXA_FILE true && ln -s $PWD/$RAW_DIR/$TAXA_FILE $TAXA_FILE )
test_file $LOCIDESC_FILE true || \
  ( test_file $RAW_DIR/$LOCIDESC_FILE true && ln -s $PWD/$RAW_DIR/$LOCIDESC_FILE $LOCIDESC_FILE )

# check scripts
test_file $SCRIPT_EVAL_LOCI true || exit_error ": script not found"
test_file $SCRIPT_BUILD_SM true || exit_error ": script not found"
test_file $SCRIPT_LOCI_CONSENSUS true || exit_error ": script not found"

echo
ask_confirm "Prefix found is \"$prefix\". Is this correct? (Y/n) " || exit_error
ask_confirm "Results directory is set to \"$RES_DIR\". Is this correct? (Y/n) " || exit_error
ask_confirm "Filters directory is set to \"$FILTERS_DIR\". Is this correct? (Y/n) " || exit_error
ask_confirm "Log directory is set to \"$LOG_DIR\". Is this correct? (Y/n) " || exit_error
ask_confirm "Supermatrices directory is set to \"$SM_DIR\". Is this correct? (Y/n) " || exit_error
ask_confirm "Per-locus consensus directory is set to \"$CONSENSUS_DIR\". Is this correct? (Y/n) " || exit_error

test_dir $LOG_DIR false || mkdir -p $LOG_DIR

if test_dir $FILTERS_DIR false; then
  filters=`ls $FILTERS_DIR/filter.* 2> /dev/null`
  n_filters=`echo $filters | wc -w`
  [ "$n_filters" == 0 ] && exit_error ": filters dir is empty"
else
  mkdir -p $FILTERS_DIR
  echo
  echo Creating filters directory in $PWD/$FILTERS_DIR
  echo Place there your custom filters with names filter.[FILTER_ID]
  echo Use $RAW_DIR/scripts/filter.cfg as template
  exit
fi

# step 1: build supermatrices

echo
cur_logfile=$LOG_DIR/build_supermatrix.log
test_dir $SM_DIR false || mkdir -p $SM_DIR
echo "building $n_filters supermatrices. Output will be appended to $cur_logfile"
for filter in $filters; do
  filtername=`echo $filter | cut -d'.' -f2-`
  echo ... filter $filtername

  mkdir -p $SM_DIR/$filtername
  #Rscript --vanilla scripts/build_supermatrix.r $LOCI_DIR $filter $prefix $SM_DIR/$filtername >> $cur_logfile
done


# step 2: evaluate all loci

  echo
  cur_logfile=$LOG_DIR/eval_loci.log
  echo "evaluating per-locus ML trees. Output will be appended to $cur_logfile"
  ask_confirm "do you want to skip already evaluated loci? (Y/n) " && force_recomp="n" || force_recomp="f"
  for i in `seq 0 $((n_kdirs-1))`; do
    echo -ne "... Block $((i+1)) of $n_kdirs\r"
    #$SCRIPT_EVAL_LOCI $i $RES_DIR $force_recomp >> $cur_logfile

    # check result files
   # for f in $LOCI_DIR/$i/$prefix.*; do
   #   f=`echo $f | rev | cut -d'/' -f1 | rev`
   #   test_file $RES_DIR/$i/$f.trees false || exit_error "Result for $i/$f is missing"
   # done
  done
  echo -ne "... Done!                      \n"

  # step 2b: build per-locus MRE consensus trees
  echo
  cur_logfile=$LOG_DIR/build_consensus_loci.log
  echo "building per-locus consensus. Output will be appended to $cur_logfile"
  mkdir -p $CONSENSUS_DIR
  for i in `seq 0 $((n_kdirs-1))`; do
    echo -ne "... Block $((i+1)) of $n_kdirs\r"
    
    output_file=$CONSENSUS_DIR/consensus.$i.trees
    test_file $output_file false && continue
    $SCRIPT_LOCI_CONSENSUS $prefix $i $RES_DIR $LOG_DIR $CONSENSUS_DIR >> $cur_logfile
  done
  echo -ne "... Done!                      \n"

  loci_consensus_file=$CONSENSUS_DIR/consensus.trees
  echo "... gathering consensus trees into $loci_consensus_file"
  rm -f $loci_consensus_file
  for i in `seq 0 $((n_kdirs-1))`; do
    test_file $CONSENSUS_DIR/consensus.$i.trees false || exit_error ": file $i not found"
    cat $CONSENSUS_DIR/consensus.$i.trees | cut -d' ' -f2 >> $loci_consensus_file
  done

  # remove support
  sed -i "s/\[[0-9]*\]//g" $loci_consensus_file

  test_file $loci_consensus_file.expanded true
  if [ "$?" == "0" ]; then
    echo "... $loci_consensus_file already expanded"
  else
    Rscript --vanilla scripts/expand_trees.r $loci_consensus_file
  fi

# step 3: evaluate supermatrices

echo
cur_logfile=$LOG_DIR/eval_supermatrices.log
echo "evaluating supermatrices. Output will be appended to $cur_logfile"
for filter in $filters; do
  filtername=`echo $filter | cut -d'.' -f2-`
  sm_dir=$SM_DIR/$filtername
  echo ... supermatrix $filtername


  # TODO: ANALYSIS HERE!

  # check result files
  test_file $sm_dir/trees.result true || exit_error
  test_file $sm_dir/trees.result.pars true || exit_error
  test_file $sm_dir/trees.result.rand true || exit_error
done

echo DONE
