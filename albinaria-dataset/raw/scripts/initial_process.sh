#!/bin/bash

# This script creates 'taxa' and %.loci.head files
# 
# We assume that taxa names in the %.loci files start with '>'
# %.loci.head file will contain the headers of %.loci but taxa names are
# converted to the indices in 'taxa' file.
#
# Run: initial_process.sh [PREFIX]

prefix=$1

# check files
loci_file=${prefix}.loci
phy_file=${prefix}.phy

#output
taxa_file=taxa
loci_head_file=${prefix}.loci.head

#generate taxa file
echo "Generate taxa file"
cut -d' ' -f 1 ${phy_file} | tail -n +2 > ${taxa_file}

n_taxa=`cat ${taxa_file} | wc -l`
echo "... there are ${n_taxa} taxa"

#generate loci head file
cut -d' ' -f 1 ${loci_file} > ${loci_head_file}

echo "Start replacing taxa"

i=1
while read taxon; do
  printf "... replace %15s ${i}/${n_taxa}\n" ${taxon}
  sed -i "s/>${taxon}$/${i}/g" ${loci_head_file}
  i=$((i+1))
done < ${taxa_file}

# check there are no unidentified taxa
test_v=`fgrep '>' ${loci_head_file} | head -n 1`
[[ ! -z ${test_v} ]] && echo "ERROR: There are unidentified taxa in loci file"
