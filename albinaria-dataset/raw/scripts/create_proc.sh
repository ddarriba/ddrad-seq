#!/bin/bash

# split .loci file into individual loci (.fas format)
#
# This script generates one msa file per locus

source scripts/aux/common.sh

# optional set of filters to evaluate
# if empty, evaluate all
shift
arg_filters=$*

# input
[ -f ${file_loci_desc} ] || { echo "Loci file ${file_loci_desc} missing"; exit; }
[ -f ${file_taxa} ] || { echo "Taxa file ${file_taxa} missing"; exit; }
[ -d ${dir_filters} ] || { echo "Filters directory ${dir_filters} missing"; exit; }

[ -f ${script_create_proc} ] || { echo "ERROR: Script ${script_stat_info} does not exist"; exit; }

# output
[ -z ${file_taxa_desc} ] && { echo "Taxa description file undefined"; exit; }
[ -z ${dir_stats} ] && { echo "Stats directory undefined"; exit; }


#-------------------------------------------------------------------------------

# search filters
if [ ! -z "${arg_filters}" ]; then
  filters=
  for filter in ${arg_filters}; do
    filters="${filters} ${dir_filters}/filter.${filter}.cfg"
  done
else
  filters=`find ${dir_filters} -name "filter.*.cfg"`
fi

n_filters=`echo $filters | wc -w`

echo "${n_filters} filters found in ${dir_filters}"

for filter in ${filters}; do
  filter_name=`echo ${filter} | rev | cut -d'.' -f2 | rev`
  echo "Processing filter ${filter_name}"

  [ -f ${filter} ] || { echo "ERROR: Filter ${filter} does not exist"; exit; }

  ARGS="${base_name} ${data_dir} ${file_loci_desc} ${dir_loci} ${filter} ${filter_name}"
  Rscript --vanilla ${script_create_proc} ${ARGS}
done

exit

