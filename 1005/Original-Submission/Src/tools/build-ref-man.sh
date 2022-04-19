#!/bin/bash
# Shell script to generate tex lines to be included in reference-manual.tex,
# one for each RMD-subroutine.
#
# Use: echo <oper> <oper>... | build-ref-man.sh <blaslevel>
#    (e.g. echo gemm symm | build-ref-man.sh 3)

echo "\pagebreak"
if [[ $1 == other ]]; then
    echo "\section{Derivatives of other subroutine(s)}"
elif [[ $1 == 4 ]]; then
    echo "\section{Adjoints of scalars}"
else
    echo "\section{Derivatives of level $1 BLAS}"
fi
read -a operation
echo "\rmdheadS{${operation[0]}}"
unset operation[0]
for op in ${operation[@]}; do
    echo "\rmdhead{$op}"
done
