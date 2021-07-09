#!/bin/bash

# common source for scripts

# directory structure
DATA_DIR=data
RESULTS_DIR=res
RAW_DATA_DIR=../raw/data

# tools
RAXML_BIN=raxmlHPC-SSE3

[ -z `which ${RAXML_BIN}` ] && { echo "ERROR: RaxML not found"; exit; }

# common arguments
subdir=$1
dir_data=${DATA_DIR}/${subdir}
dir_raw_data=${RAW_DATA_DIR}/${subdir}
dir_loci=${dir_data}/loci

dir_results=${RESULTS_DIR}/${subdir}

# ignore used arguments for scripts
shift

[ -z "$subdir" -o -z "$dir_results" ] && { echo "ERROR: Script requires subdir + results_dir arguments"; exit; }
[ -d $dir_data ] || { echo "ERROR: Data directory $dir_data does not exist"; exit; }
[ -d $RESULTS_DIR ] ||  { echo "ERROR: Results directory $RESULTS_DIR does not exist"; exit; }
[ -d $dir_loci ] || { echo "ERROR: Loci directory $dir_loci does not exist"; exit; }
[ -d $dir_raw_data ] || { echo "ERROR: Raw data directory $dir_raw_data does not exist"; exit; }

mkdir -p ${dir_results}

# useful files
file_loci_desc=`find ${dir_raw_data} -name "*.loci.desc"`
base_name=`echo $file_loci_desc | rev | cut -d'/' -f1 | cut -d'.' -f3- | rev`
file_seeds=scripts/seeds

[ -z $file_loci_desc ] && { echo "ERROR: Loci desc file in $file_loci_desc not found"; exit; }
[ -f $file_seeds ] || { echo "ERROR: Seeds file $file_seeds not found"; exit; }

file_taxa=${dir_raw_data}/taxa
file_taxa_desc=${dir_raw_data}/taxa.desc
file_loci_head=${dir_raw_data}/${base_name}.loci.head
file_loci_desc=${dir_raw_data}/${base_name}.loci.desc

file_tmp=${data_dir}/tempfile

dir_config=${dir_raw_data}/config
dir_stats=${dir_raw_data}/stats
dir_filters=${dir_config}/filters

[ -d $dir_config ] || { echo "ERROR: Config directory $dir_config does not exist"; exit; }
[ -d $dir_filters ] || { echo "ERROR: Filters directory $dir_filters does not exist"; exit; }

echo "--- Preprocessing info ---"
echo "Data directory:    $dir_data"
echo "Results directory: $dir_results"
echo "Base name:         $base_name"

# scripts
script_stat_info=scripts/aux/stat_info.r
script_create_proc=scripts/aux/create_proc.r
script_build_supermatrix=scripts/aux/build_supermatrix.r
