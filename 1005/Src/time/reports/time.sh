#!/bin/bash
# Usage examples, assuming $HOST is tg
#    echo tg openblas 4 100 | time.sh      (timermd -t4 -n100)
#    echo . mkl 4 1000 10 100M | time.sh   (. matches all hosts)
#    time.sh < timepar.dat                 (only run lines beginning with tg)
#    time.sh << EOF
#      -s2 -m100M                          (add -s2 -m100M to timermd)
#      tg    mkl 1 100                     (current host is tg so run this line)
#      tg    mkl 4 100                     (and this)
#      mimir mkl 1 100                     (but not this)
#      mimir mkl 4 100                     (or this)
#    EOF
#    time.sh -f < timepar.dat              (fast: don't make, don't add flags)

while read line; do
    if [[ $line =~ ^- && $1 != -f ]]
    then flags=$line
    else read computer blas t n s m <<< $line
    fi
    
    if [[ $computer = $HOST || $computer = . ]]; then
        echo
        echo $HOST:$blas
        if [[ $1 != -f ]]; then 
            echo Making...
            cd ~/blas_rmd
            make clean > make.log
            make -j time blas=$blas optimize=1 >> make.log
        fi
            echo Running timermd -t$t -n$n $flags
        echo ...and copying results to reports/$blas-$HOST.dat
        cd ~/blas_rmd/time
        timermd -t$t -n$n $flags | tee reports/$blas-$HOST.dat
    fi
done
