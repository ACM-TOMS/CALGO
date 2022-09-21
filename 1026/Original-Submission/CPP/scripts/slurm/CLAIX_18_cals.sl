#!/usr/local_rwth/bin/zsh

#SBATCH --job-name=CALS-EXP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=40G
#SBATCH --output=output.%J.txt
#SBATCH --time=05:00:00
#SBATCH --partition=c18m
#SBATCH --exclusive

# if NOT CTF
#module load clang/9.0
# else
module load gcc/9
#module switch intelmpi/2018 openmpi/4.0.3
# fi

#module load LIBRARIES
#module load intelmkl/2020

source ~/.zshrc.local
cd ${CALS_DIR} || exit

if [ -d "build" ]; then
  rm -rvf build/*
  rm -rvf build
  mkdir build
  touch build/.gitkeep
else
  mkdir build
fi

cd build || exit
#CC=clang CXX=clang++ ./../extern/cmake/bin/cmake -DCMAKE_BUILD_TYPE=Release -DWITH_MKL=On -DWITH_EXPERIMENTS=Off -DWITH_DIAGNOSTICS=1 ..
CC=gcc CXX=g++ ./../extern/cmake/bin/cmake -DCMAKE_BUILD_TYPE=Release -DWITH_MKL=On -DWITH_EXPERIMENTS=Off -DWITH_DIAGNOSTICS=1 ..
make -j 48

numactl -H
lscpu
#export GOMP_CPU_AFFINITY="0 1 2 6 7 8 12 13 14 18 19 20 3 4 5 9 10 11 15 16 17 21 22 23"  # CLAIX18 xeon platinum 8160

numactl --cpubind=0,1 --membind=0,1 -- numactl -show
export OMP_NUM_THREADS=24
export OMP_THREAD_LIMIT=24
numactl --cpubind=0,1 --membind=0,1 -- ./src/experiments/experiments_MKL 24
#numactl --cpubind=0,1 --membind=0,1 -- ./src/experiments/benchmark_cals_mttkrp_MKL 24
#numactl --cpubind=0,1 --membind=0,1 -- ./src/experiments/benchmark_other_mttkrp_MKL 24

#numactl --cpubind=0 --membind=0 -- numactl -show
#export OMP_NUM_THREADS=12
#export OMP_THREAD_LIMIT=12
#numactl --cpubind=0 --membind=0 -- ./src/experiments/experiments_MKL 12
#numactl --cpubind=0 --membind=0 -- ./src/experiments/benchmark_cals_mttkrp_MKL 12
#numactl --cpubind=0 --membind=0 -- ./src/experiments/benchmark_other_mttkrp_MKL 12


numactl --cpubind=0 --membind=0 -- numactl -show
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1
numactl --cpubind=0 --membind=0 -- ./src/experiments/experiments_MKL 1
#numactl --cpubind=0 --membind=0 -- ./src/experiments/benchmark_cals_mttkrp_MKL 1
#numactl --cpubind=0 --membind=0 -- ./src/experiments/benchmark_other_mttkrp_MKL 1

# make clean
