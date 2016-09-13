base_name=$1
loci_file=${base_name}.loci.head

list=
locus_id=1
while read l; do
  if [[ $l == "//" ]]; then
    printf "%8s " $locus_id
    for i in `seq 1 98`; do
      [[ $list =~ (^| )${i}($| ) ]] && printf 1 || printf 0
    done
    printf '\n'
    list=
    locus_id=$((locus_id + 1))
  else
    list+=" $l"
  fi
done < $loci_file
