#!/bin/bash

# Usage:
#
#  eval_loci.sh SUBDIR K_DIR
#
#    SUBDIR = data set subdirectory in data/
#    K_DIR = 0..(n_kdirs-1)
#

# static configuration

n_parsimony=5
n_random=5
n_replicates=$((n_parsimony + n_random))

error_folder=error

source scripts/aux/common.sh

function move_results {
  # expect $1: execution name
	# expect $2: output directory
  #echo "Appending $1 to $2"
  cat RAxML_bestTree.$1 >> $2
  rm RAxML_*.$1
}

function move_results {
  # expect $1: execution name
  mkdir -p ${error_folder}
  mv RAxML_bestTree.$1 ${error_folder}
}

function process {
  # expect $1: raxml
	# expect $2: dataset
	# expect $3: result base path
	# expect $4: # replicates from parsimony tree
  # expect $5: # replicates from random tree

  threads_per_process=1
  raxml_cmd=$1
  in_path=$2
  out_path=$3
  n_parsimony_local=$4
  n_random_local=$5
  log_path=$6

  if [[ $threads_per_process == "1" ]]; then
    cmd_args=
  else
    cmd_args="-T $threads_per_process"
  fi

  # clean output
  [ -f $out_path ] && rm $out_path

  counter=1
  while read seed; do
    exec_name=`echo $in_path | rev | cut -d'.' -f1 | rev`
    exec_name+=.$seed
		cmd="$raxml_cmd -s $in_path -n $exec_name -p $seed -m GTRGAMMA $cmd_args"
    if [[ $n_parsimony_local > 0 ]]; then
      n_parsimony_local=$((n_parsimony_local - 1))
    else
      cmd="$cmd -d"
      n_random_local=$((n_random_local - 1))
      [ $n_random_local -lt 0 ] && break
    fi

	  echo -ne "$counter [$seed]        \r"

    # clear output (if exists)
    [ -f RAxML_info.${exec_name} ] && rm RAxML_*.${exec_name}
    eval $cmd >> $log_path && move_results $exec_name $out_path || move_error $exec_name

    counter=$((counter+1))
	done < $file_seeds
}

loci_subset=$1
[ -z ${loci_subset} ] && { echo "Loci subset missing"; exit; }

	datFolder=${dir_loci}/${loci_subset}
	resFolder=${dir_results}/${loci_subset}
  force_recomp=$3
	logFolder=log/$1
  mkdir -p $resFolder || exit
  mkdir -p $logFolder || exit

  n_loci=`ls ${datFolder}/${base_name}.locus.* | wc -l`

  echo "Found $n_loci samples in $datFolder"

  datFN=`echo $datFolder | rev | cut -d'/' -f1 | rev`

  echo "$seedlist" > $resFolder/${loci_subset}.seeds
	for datFile in $datFolder/*; do
		datName=$(basename "$datFile")
    file_trees="$resFolder/$datName.trees"
    [ -f $file_trees ] && n_trees=`wc -l $file_trees | cut -d' ' -f1` || n_trees=0
    if [ "$force_recomp" == "f" -o ${n_trees} != ${n_replicates} ]; then
  		echo `date` "Processing $datFile -> $resFolder/$datName";
  		process $RAXML_BIN $datFile "$resFolder/$datName.trees" $n_parsimony $n_random "$logFolder/$datName.log";
    fi
	done;

