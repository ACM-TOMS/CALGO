#!/usr/local_rwth/bin/zsh

#SBATCH --job-name=Bench-CU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=output.%J.txt
#SBATCH --time=02:00:00
#SBATCH --partition=c18g
#SBATCH --gres=gpu:volta:2
#SBATCH --exclusive
# comment -B 1:1

module load DEVELOP
module load gcc/9
module load cuda/110

source ~/.zshrc.local
cd ${CALS_DIR} || exit

if [ -d "build_cuda" ]; then
  rm -rvf build_cuda/*
  rm -rv build_cuda
  mkdir build_cuda
else
  mkdir build_cuda
fi

cd build_cuda || exit
#CC=clang CXX=clang++ ./../extern/cmake/bin/cmake -DCMAKE_BUILD_TYPE=Release -DWITH_MKL=On -DWITH_EXPERIMENTS=Off -DWITH_DIAGNOSTICS=2 -DWITH_CUBLAS=On ..
CC=gcc CXX=g++ ./../extern/cmake/bin/cmake -DCMAKE_BUILD_TYPE=Release -DWITH_MKL=On -DWITH_EXPERIMENTS=Off -DWITH_DIAGNOSTICS=2 -DWITH_CUBLAS=On ..
make -j 48

numactl -H
numactl --cpubind=0,1 --membind=0,1 -- numactl -show
export OMP_NUM_THREADS=24
export OMP_THREAD_LIMIT=24
numactl --cpubind=0,1 --membind=0,1 -- ./src/experiments/experiments_MKL 24
#numactl --cpubind=0 --membind=0 -- ./src/experiments/benchmark_cals_mttkrp_MKL 24

# make clean
