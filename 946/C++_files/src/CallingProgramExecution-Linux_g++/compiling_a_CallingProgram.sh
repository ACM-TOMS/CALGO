#!/bin/sh


#TO CHANGE #################################
reliadiffdir=$PWD/..
############################################

echo RELIADIFF directory : $reliadiffdir

echo '############################################'
echo Cleaning previous built versions...
echo '############################################'
echo
make clean -C $reliadiffdir
echo
echo '############################################'
echo Creating RELIADIFF library...
echo '############################################'
echo
make -C $reliadiffdir
echo

echo '########################################################'
echo Compiling, linking and executing the Calling Program...
echo '########################################################'
echo
g++ -ansi -Wall -pedantic -g -I$reliadiffdir -c ../SimpleCallingProgramExample.c 
g++ -ansi -Wall -pedantic -g -o calling_program SimpleCallingProgramExample.o -lm  -L$reliadiffdir -lreliadiff
rm SimpleCallingProgramExample.o

echo '############################################'
echo Output of the Calling Program
echo '############################################'
echo

./calling_program