#!/bin/bash

# This script creates 'taxa' and %.loci.head files
# 
# We assume that taxa names in the %.loci files start with '>'
# %.loci.head file will contain the headers of %.loci but taxa names are
# converted to the indices in 'taxa' file.
#
# Run: initial_process.sh [DATA_DIR]

source scripts/aux/common.sh

# input
[[ -f ${file_loci} ]] || { echo "Loci file missing"; exit; }
[[ -f ${file_phylip} ]] || { echo "PHYLIP file missing"; exit; }

# output
[[ -z ${file_taxa} ]] && { echo "Taxa file undefined"; exit; }
[[ -z ${file_loci_head} ]] && { echo "Loci head file undefined"; exit; }

#generate taxa file
echo "Generate taxa file"
cut -d' ' -f 1 ${file_phylip} | tail -n +2 > ${file_taxa}

n_taxa=`cat ${file_taxa} | wc -l`
echo "... there are ${n_taxa} taxa"

#generate loci head file
cut -d' ' -f 1 ${file_loci} > ${file_loci_head}

echo "Start replacing taxa"

i=1
while read taxon; do
  printf "... replace %15s ${i}/${n_taxa}\n" ${taxon}
  sed -i "s/>${taxon}$/${i}/g" ${file_loci_head}
  i=$((i+1))
done < ${file_taxa}

# check there are no unidentified taxa
test_v=`fgrep '>' ${file_loci_head} | head -n 1`
[[ ! -z ${test_v} ]] && echo "ERROR: There are unidentified taxa in loci file"
