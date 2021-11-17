while read computer; do
    scp $computer:blas_rmd/time/reports/\*-$computer.dat .
done <<EOF
mimir
jotunn
ubuntu
EOF
