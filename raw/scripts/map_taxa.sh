#!/bin/bash

source scripts/aux/common.sh

#input
[[ -f ${file_loci_head} ]] || { echo "Loci file missing"; exit; }

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
done < ${file_loci_head}
