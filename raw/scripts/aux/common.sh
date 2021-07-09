#!/bin/bash

# common source for scripts

# directory structure
DATA_DIR=data
FILES_SUBDIR=outfiles

# tools
RAXML_BIN=raxmlHPC-SSE3

[ -z `which ${RAXML_BIN}` ] && { echo "ERROR: RaxML not found"; exit; }

# common arguments
subdir=$1
data_dir=${DATA_DIR}/$1
input_dir=${data_dir}/${FILES_SUBDIR}

[ -z $subdir ] && { echo "ERROR: Script requires 1 argument"; exit; }
[ -d $data_dir ] || { echo "ERROR: Data directory $data_dir does not exist"; exit; }
[ -d $input_dir ] || { echo "ERROR: Input data directory $input_dir does not exist"; exit; }

# useful files
phy_file=`find $input_dir -name "*.phy"`
base_name=`echo $phy_file | rev | cut -d'/' -f1 | cut -d'.' -f2- | rev`

[ -z $phy_file ] && { echo "ERROR: PHYLIP file in $input_dir not found"; exit; }

echo "--- Preprocessing info ---"
echo "Data directory: $data_dir"
echo "Base name:      $base_name"

file_loci=${input_dir}/${base_name}.loci
file_phylip=${input_dir}/${base_name}.phy
file_taxa=${data_dir}/taxa
file_taxa_desc=${data_dir}/taxa.desc
file_loci_head=${data_dir}/${base_name}.loci.head
file_loci_desc=${data_dir}/${base_name}.loci.desc

file_tmp=${data_dir}/tempfile

dir_loci=${data_dir}/loci
dir_config=${data_dir}/config
dir_stats=${data_dir}/stats
dir_filters=${dir_config}/filters

[ -d $dir_config ] || { echo "ERROR: Config directory $dir_config does not exist"; exit; }
[ -d $dir_filters ] || { echo "ERROR: Filters directory $dir_filters does not exist"; exit; }


# scripts
script_stat_info=scripts/aux/stat_info.r
script_create_proc=scripts/aux/create_proc.r
script_build_supermatrix=scripts/aux/build_supermatrix.r
