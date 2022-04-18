#!/bin/bash
# run with:
#   make clean
#   ./timedemo2.sh
#   mv demo2out.txt ../time/reports
set -x
export slicot=1 BLAS=mkl optimize=1
cd ..
#make clean
make -j demo
cd demo
rv=(100 100 500)
nv=(1000 5000 1000)
> demo2out.txt
for T in 1 8; do
    export MKL_NUM_THREADS=$T
    echo $T threads: >> demo2out.txt
    for i in 0 1 2; do
        r=${rv[$i]} n=${nv[$i]}
        { time -p demo2_driver -n $r $n; } 2> >(tee -a demo2out.txt) 
        { time -p demo2_driver $r $n; } 2> >(tee -a demo2out.txt)
    done
done
