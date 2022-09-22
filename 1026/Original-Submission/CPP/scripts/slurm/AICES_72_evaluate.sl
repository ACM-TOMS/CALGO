#!/usr/local_rwth/bin/zsh

#SBATCH --job-name=EVAL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err
#SBATCH --time=00:10:00
#SBATCH -p ih
#SBATCH -A aices2
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
  rm -rvf build/*
  rm -rvf build
  mkdir build
else
  mkdir build
fi

cd build || exit
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_TESTS=Off -DWITH_MKL=On -DWITH_OPENBLAS=Off -DWITH_BLIS=Off -DWITH_CUBLAS=Off ..
make -j 24 Evaluator_MKL
#               Threads  MHz(AVX) (D)FPC
./Evaluator_MKL       1      2500     16
./Evaluator_MKL       12     2500     16

make clean
