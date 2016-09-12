#!/bin/sh
# split .loci file into individual loci (.fas format)

base_name=Albinaria98inds_c91d6m4p3
loci_file=${base_name}.loci
output_dir=perlocus
mkdir -p $output_dir

tmp_file=$output_dir/tmp

loci_sections=`fgrep -n "//" $loci_file | cut -d':' -f 1`
cur_line=1
locus_id=1
n_loci=`echo $loci_sections | wc -w`
for locus_end in $loci_sections; do
  subdir=${output_dir}/$((locus_id / 1000))
  mkdir -p $subdir
  out_file=${subdir}/${base_name}.locus.${locus_id}
  echo $locus_id / $n_loci
  locus_id=$((locus_id+1))

  loop_start=$cur_line
  cur_line=$((locus_end+1))
  # skip existing files
  [ -f $out_file ] && continue

  # use tmp_file for assuring integrity of out_file
  sed -n "${loop_start},$((locus_end-1))p" $loci_file > $tmp_file
  sed -i "s# \+#\n#g" $tmp_file
  mv $tmp_file $out_file
done
