#!/bin/bash
for file in *.dat; do
    IFS=".-" read blas computer rest <<< "$file"
    awk -v blas=$blas -v computer=$computer  '
        /^The table/ {pon=0}
        pon          {print computer,blas,$3,$4,$1,$5,$6,$9,$10} 
        /Operation/  {pon=1}' $file
done \
    | sed 's/openblas/open/' \
    | sed 's/accelerate/accel/' \
    | sort -k1,2 -k3,3n -k4,4n > times.txt
