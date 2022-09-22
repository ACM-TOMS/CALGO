#!/bin/bash

#SBATCH -J CALS
#SBATCH -t 48:00:00
#SBATCH -p ih
#SBATCH -A aices2
#SBATCH --mem=58G

#SBATCH -o cals.log
#SBATCH -C hwx2680
#SBATCH --nodelist=linuxihdc072
#SBATCH --exclusive

module purge
module load DEVELOP
module load gcc/8
module load LIBRARIES
module load intelmkl/2019
module load cmake/3.13.2

source ~/.zshrc.local
cd ${CALS_DIR} || exit

export GOMP_CPU_AFFINITY="0,2,4,6,8,10,12,14,16,18,20,22"  # Node 72 has even number core on socket 0 (12 total)

if [ -d "build" ]; then
  rm -rf build/*
  rmdir build
  mkdir build
else
  mkdir build
fi

cd build || exit
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_TESTS=Off -DWITH_MKL=On -DWITH_OPENBLAS=Off -DWITH_BLIS=Off -DWITH_CUBLAS=Off ..
make -j 24
#./Experiment_MKL 12
#./Experiment_MKL 1

#./Experiment_OPENBLAS 1
#./Experiment_OPENBLAS 12

#./Experiment_BLIS 1
#./Experiment_BLIS 12
make clean
