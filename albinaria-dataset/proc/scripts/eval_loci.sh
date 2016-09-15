#!/bin/bash

raxml_bin="raxmlHPC-SSE3"
seedlist_file=scripts/seeds
n_parsimony=50
n_random=50
n_replicates=$((n_parsimony + n_random))

function move_results {
  # expect $1: execution name
	# expect $2: output directory
  #echo "Appending $1 to $2"
  cat RAxML_bestTree.$1 >> $2
  #echo Removing RAxML_*.$1
  rm RAxML_*.$1
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
  n_parsimony=$4
  n_random=$5
  log_path=$6

  if [[ $threads_per_process == "1" ]]; then
    cmd_args=
  else
    cmd_args="-T $threads_per_process"
  fi

  # clean output
  if [ -f $out_path ]; then
    rm $out_path
  fi

  while read seed; do
    exec_name=`echo $in_path | rev | cut -d'.' -f1 | rev`
    exec_name+=.$seed
		cmd="$raxml_cmd -s $in_path -n $exec_name -p $seed -m GTRGAMMA $cmd_args"
    if [[ $n_parsimony > 0 ]]; then
      n_parsimony=$((n_parsimony - 1))
    else
      cmd="$cmd -d"
      n_random=$((n_random - 1))
    fi
		#echo "$cmd"
    eval $cmd >> $log_path && move_results $exec_name $out_path;
	done < $seedlist_file
}

[[ $1 == "" ]] && echo "Input dir missing" && exit
[[ $2 == "" ]] && echo "Results dir missing" && exit

test_bin=`which $raxml_bin`
if [[ -f $test_bin ]]; then
	datFolder=loci/$1
	resFolder=$2/$1
	logFolder=log/$1
  mkdir -p $resFolder || exit
  mkdir -p $logFolder || exit
  datFN=`echo $datFolder | rev | cut -d'/' -f1 | rev`
  echo "$seedlist" > $resFolder/$datFN.seeds
	for datFile in $datFolder/*; do
		datName=$(basename "$datFile")
		echo `date` "Processing $datFile -> $resFolder/$datName";
		process $raxml_bin $datFile "$resFolder/$datName.trees" $n_parsimony $n_random "$logFolder/$datName.log";
	done;
else
	echo "Error: Binary $raxml_bin does not exist."
fi

