#!/usr/local_rwth/bin/zsh

#SBATCH --job-name=CALS_MAT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --output=cals-matlab-%j.out
#SBATCH --error=cals-matlab-%j.err
#SBATCH --time=8:00:00
#SBATCH --partition=c18m
#SBATCH --exclusive

module purge
module load DEVELOP
module load MISC
module load matlab/2019b

source ~/.zshrc.local
#cd ${CALS_DIR} || exit
#
#if [ -d "build_matlab" ]; then
#  rm -rvf build_matlab/*
#  rm -rvf build_matlab
#  mkdir build_matlab
#else
#  mkdir build_matlab
#fi
#
#cd build_matlab || exit
#cmake -DCMAKE_BUILD_TYPE=Release -DWITH_MATLAB=On ..
#make -j 48

exppath="${CALS_DIR}/matlab/matlab_src"

#export GOMP_CPU_AFFINITY="0 1 2 6 7 8 12 13 14 18 19 20 3 4 5 9 10 11 15 16 17 21 22 23"  # CLAIX18 xeon platinum 8160
export MATLAB_LOG_DIR=.
export MATLABPATH=$exppath

cd $exppath || exit
pwd

export MKL_NUM_THREADS=24
export OMP_NUM_THREADS=24
matlab -nodisplay -nodesktop -nosplash -logfile /dev/null <<EOF
addpath('${exppath}');
TTB_experiment;
quit();
EOF

#export MKL_NUM_THREADS=12
#export OMP_NUM_THREADS=12
#matlab -nodisplay -nodesktop -nosplash -logfile /dev/null <<EOF
#addpath('${exppath}');
#TTB_benchmark;
#quit();
#EOF

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
matlab -singleCompThread -nodisplay -nodesktop -nosplash -logfile /dev/null <<EOF
addpath('${exppath}');
TTB_experiment;
quit();
EOF
