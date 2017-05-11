#!/bin/bash

# split .loci file into individual loci (.fas format)
#
# This script generates one msa file per locus

source scripts/aux/common.sh

# input
[ -f ${file_loci} ] || { echo "Loci file missing"; exit; }

# output
[ -z ${dir_loci} ] && { echo "Loci directory undefined"; exit; }

output_dir=$data_dir/loci
mkdir -p ${dir_loci}

loci_sections=`fgrep -n "//" $file_loci | cut -d':' -f 1`
cur_line=1
locus_id=1
n_loci=`echo $loci_sections | wc -w`
for locus_end in $loci_sections; do
  subdir=${dir_loci}/$((locus_id / 1000))
  mkdir -p $subdir
  out_file=${subdir}/${base_name}.locus.${locus_id}
  echo $locus_id / $n_loci
  locus_id=$((locus_id+1))

  loop_start=$cur_line
  cur_line=$((locus_end+1))

  # skip existing files
  [ -f $out_file ] && continue

  # use a temporary for assuring integrity of out_file
  sed -n "${loop_start},$((locus_end-1))p" $file_loci > ${file_tmp}
  sed -i "s# \+#\n#g" ${file_tmp}
  mv ${file_tmp} $out_file
done
