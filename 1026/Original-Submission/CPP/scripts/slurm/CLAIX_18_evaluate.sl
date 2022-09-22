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
#SBATCH --partition=c18m
#SBATCH --account=rwth0575
#SBATCH --exclusive

module purge
module load DEVELOP
module load gcc/8
module load LIBRARIES
module load intelmkl/2019
module load cmake/3.13.2

source ~/.zshrc.local
cd ${CALS_DIR} || exit
export GOMP_CPU_AFFINITY="0 1 2 6 7 8 12 13 14 18 19 20 3 4 5 9 10 11 15 16 17 21 22 23"  # CLAIX18 xeon platinum 8160

if [ -d "build" ]; then
  rm -rvf build/*
  rm -rvf build
  mkdir build
  touch build/.gitkeep
else
  mkdir build
fi

cd build || exit
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_TESTS=Off -DWITH_MKL=On -DWITH_OPENBLAS=Off -DWITH_BLIS=Off -DWITH_CUBLAS=Off ..
make -j 48 Evaluator_MKL
#               Threads  MHz(AVX) (D)FPC Size
./Evaluator_MKL       1      3500     32 1000
./Evaluator_MKL       12     2600     32 1000
./Evaluator_MKL       24     2000     32 1000

./Evaluator_MKL       1      3500     32 1921
./Evaluator_MKL       12     2600     32 1921
./Evaluator_MKL       24     2000     32 1921

./Evaluator_MKL       1      3500     32 2828
./Evaluator_MKL       12     2600     32 2828
./Evaluator_MKL       24     2000     32 2828

./Evaluator_MKL       1      3500     32 5196
./Evaluator_MKL       12     2600     32 5196
./Evaluator_MKL       24     2000     32 5196

make clean
