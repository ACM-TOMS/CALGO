#!/bin/zsh

source ~/.zshrc.local
cd ${CALS_DIR} || exit
export GOMP_CPU_AFFINITY="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19"  # CLAIX18 xeon platinum 8160

if [ -d "build" ]; then
  rm -rvf build/*
  rm -rvf build
  mkdir build
  touch build/.gitkeep
else
  mkdir build
fi

cd build || exit
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_MKL=On ..
make -j 40 .
./experiment_jackkniffing_MKL 10
