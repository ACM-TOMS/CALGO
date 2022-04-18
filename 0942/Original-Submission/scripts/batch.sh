#!/bin/bash

PREFIX=`cat $1 | grep PREFIX | cut -d " " -f 3-`
BIN=`cat $1 | grep BIN_PATH | awk '{print $3}'`
RES=`cat $1 | grep RESULS_PATH | awk '{print $3}'`
MACHINE=`cat $1 | grep MACHINE_NAME | awk '{print $3}'`
THREADS=`cat $1 | grep OMP_THREADS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`


# Naive
LENGTHS=`cat $1 | grep naive_LENGTHS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
TIMESTEPS=`cat $1 | grep naive_TIMESTEPS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
SIZES=`cat $1 | grep naive_SIZES | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
BENCHS=`cat $1 | grep naive_BENCHS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`

for BENCH in $BENCHS ; do
  for LENGTH in $LENGTHS ; do
    for TIMESTEP in $TIMESTEPS ; do
      FILE="$RES/$MACHINE.$BENCH.t$TIMESTEP.l$LENGTH"
      for THREAD in $THREADS ; do
        echo "OMP_NUM_THREADS=$THREAD" 1>> $FILE
        for SIZE in $SIZES ; do
          echo "OMP_NUM_THREADS=$THREAD $PREFIX $BENCH $SIZE $SIZE $SIZE $SIZE $SIZE $SIZE $TIMESTEP $LENGTH"
          OMP_NUM_THREADS=$THREAD $PREFIX $BIN/$BENCH $SIZE $SIZE $SIZE $SIZE $SIZE $SIZE $TIMESTEP $LENGTH 2> /dev/null 1>> $FILE
        done
        echo " " 1>> $FILE
      done
    done
  done
done


# Rivera
LENGTHS=`cat $1 | grep rivera_LENGTHS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
TIMESTEPS=`cat $1 | grep rivera_TIMESTEPS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
SIZES=`cat $1 | grep rivera_SIZES | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
BENCHS=`cat $1 | grep rivera_BENCHS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
BLOCKS=`cat $1 | grep rivera_BLOCKS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`

for BENCH in $BENCHS ; do
  for LENGTH in $LENGTHS ; do
    for TIMESTEP in $TIMESTEPS ; do
      FILE="$RES/$MACHINE.$BENCH.t$TIMESTEP.l$LENGTH"
      for THREAD in $THREADS ; do
        echo "OMP_NUM_THREADS=$THREAD" 1>> $FILE
        for BLOCKJ in $BLOCKS ; do
#          for BLOCKI in $BLOCKS ; do
            for SIZE in $SIZES ; do
              echo "OMP_NUM_THREADS=$THREAD $PREFIX $BENCH $SIZE $SIZE $SIZE $SIZE $BLOCKJ $SIZE $TIMESTEP $LENGTH"
              OMP_NUM_THREADS=$THREAD $PREFIX $BIN/$BENCH $SIZE $SIZE $SIZE $SIZE $BLOCKJ $SIZE $TIMESTEP $LENGTH 2> /dev/null 1>> $FILE
            done
#          done
        done
        echo " " 1>> $FILE
      done
    done
  done
done


# Time-skewing
LENGTHS=`cat $1 | grep time_LENGTHS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
TIMESTEPS=`cat $1 | grep time_TIMESTEPS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
SIZES=`cat $1 | grep time_SIZES | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
BENCHS=`cat $1 | grep time_BENCHS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
BLOCKS=`cat $1 | grep time_BLOCKS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`

for BENCH in $BENCHS ; do
  for LENGTH in $LENGTHS ; do
    for TIMESTEP in $TIMESTEPS ; do
      FILE="$RES/$MACHINE.$BENCH.t$TIMESTEP.l$LENGTH"
      for THREAD in $THREADS ; do
        echo "OMP_NUM_THREADS=$THREAD" 1>> $FILE
        for BLOCKJ in $BLOCKS ; do
#          for BLOCKI in $BLOCKS ; do
            for SIZE in $SIZES ; do
              echo "OMP_NUM_THREADS=$THREAD $PREFIX $BENCH $SIZE $SIZE $SIZE $SIZE $BLOCKJ $SIZE $TIMESTEP $LENGTH"
              OMP_NUM_THREADS=$THREAD $PREFIX $BIN/$BENCH $SIZE $SIZE $SIZE $SIZE $BLOCKJ $SIZE $TIMESTEP $LENGTH 2> /dev/null 1>> $FILE
            done
#          done
        done
        echo " " 1>> $FILE
      done
    done
  done
done


# Oblivious
LENGTHS=`cat $1 | grep oblivious_LENGTHS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
TIMESTEPS=`cat $1 | grep oblivious_TIMESTEPS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
SIZES=`cat $1 | grep oblivious_SIZES | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
BENCHS=`cat $1 | grep oblivious_BENCHS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`
BLOCKS=`cat $1 | grep oblivious_BLOCKS | awk '{ for(i=3;i<=NF;i++) {print $i} }'`

for BENCH in $BENCHS ; do
  for LENGTH in $LENGTHS ; do
    for TIMESTEP in $TIMESTEPS ; do
      FILE="$RES/$MACHINE.$BENCH.t$TIMESTEP.l$LENGTH"
      for THREAD in $THREADS ; do
        echo "OMP_NUM_THREADS=$THREAD" 1>> $FILE
        for BLOCKI in $BLOCKS ; do
          for SIZE in $SIZES ; do
            echo "OMP_NUM_THREADS=$THREAD $PREFIX $BENCH $SIZE $SIZE $SIZE $BLOCKI $SIZE $SIZE $TIMESTEP $LENGTH"
            OMP_NUM_THREADS=$THREAD $PREFIX $BIN/$BENCH $SIZE $SIZE $SIZE $BLOCKI $SIZE $SIZE $TIMESTEP $LENGTH 2> /dev/null 1>> $FILE
          done
        done
        echo " " 1>> $FILE
      done
    done
  done
done

