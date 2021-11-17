#!/bin/bash

cmakeExists=`which cmake | wc -l`
#echo $cmakeExists
if [ "$cmakeExists" = "0" ]
then
    echo "[31m[1mWarning: to correctly build AGSol, the 'cmake' package must be installed.[0m"
    echo
    exit
fi

gplusplusExists=`which g++ | wc -l`
#echo $gplusplusExists
if [ "$gplusplusExists" = "0" ]
then
    echo "[31m[1mWarning: to correctly build AGSol, the 'g++' package must be installed.[0m"
    echo
    exit
fi

gmpExists=`ldconfig -p | grep libgmp | wc -l`
#echo $gmpExists
if [ "$gmpExists" = "0" ]
then
    echo "[31m[1mWarning: to correctly build AGSol, the 'libgmp' package must be installed.[0m"
    echo
    exit
fi

mpfrExists=`ldconfig -p | grep libmpfr | wc -l`
#echo $mpfrExists
if [ "$mpfrExists" = "0" ]
then
    echo "[31m[1mWarning: to correctly build AGSol, the 'libmpfr' package must be installed.[0m"
    echo
    exit
fi

boostExists=`ldconfig -p | grep libboost | wc -l`
#echo $boostExists
if [ "$boostExists" = "0" ]
then
    echo "[31m[1mWarning: to correctly build AGSol, the 'libboost' package must be installed.[0m"
    echo
    exit
fi

cgalExists=`ldconfig -p | grep libCGAL | wc -l`
#echo $cgalExists
if [ "$cgalExists" = "0" ]
then
    echo "[31m[1mWarning: to correctly build AGSol, the 'libCGAL' package must be installed.[0m"
    echo
    exit
fi

while test $# -gt 0
do
    useSolver=`echo $1 | grep 'DXPRESS_PATH=\|DGLPK_PATH=' | wc -l`

    if [ "$useSolver" -eq "1" ]
    then
        param=$param" "$1

        isXpress=`echo $1 | grep 'DXPRESS_PATH=' | wc -l`
    
        if [ $isXpress = "1" ]
        then
            path=`echo $1 | sed 's/\(.*\)\=\(.*\)/\2/'`
           
            if [ ! -d $path ]
            then
                echo "[31m[1mWarning: the path '"$path"' does not exists. Please verify the correct place where Xpress is installed.[0m"
                echo
                exit
            fi
            
            if [ ! -d $path/lib ]
            then
                echo "[31m[1mWarning: the path '"$path"/lib' does not exists. Please verify the correct place where Xpress is installed.[0m"
                echo
                exit
            fi

            xpressExists=`ls $path/lib | grep libxprs | wc -l`
            #echo $xpressExists 
            if [ "$xpressExists" = "0" ]
            then
                echo "[31m[1mWarning: there is no Xpress library in the path '"$path"/lib'. Please verify the correct place where Xpress is installed.[0m"
                echo
                exit
            fi
        else
            path=`echo $1 | sed 's/\(.*\)\=\(.*\)/\2/'`
           
            if [ ! -d $path ]
            then
                echo "[31m[1mWarning: the path '"$path"' does not exists. Please verify the correct place where Glpk is installed.[0m"
                echo
                exit
            fi
            
            if [ ! -d $path/lib ]
            then
                echo "[31m[1mWarning: the path '"$path"/lib' does not exists. Please verify the correct place where Glpk is installed.[0m"
                echo
                exit
            fi

            glpkExists=`ls $path/lib | grep libglpk | wc -l`
            #echo $glpkExists 
            if [ "$glpkExists" = "0" ]
            then
                echo "[31m[1mWarning: there is no Gplk library in the path '"$path"/lib'. Please verify the correct place where Glpk is installed.[0m"
                echo
                exit
            fi
        fi
    fi

    shift
done

if [ "$param" == "" ]
then
    echo "[31m[1m?> ./build-all.sh -DGLPK_PATH=<glpk path*> -DXPRESS_PATH=<xpress path**>[0m"
    #echo "[31m[1m* glpk library must be in '<glpk path>/lib' and glpk headers must be in '<glpk path>/include'[0m"
    #echo "[31m[1m** xpress library must be in '<xpress path>/lib' and xpress headers must be in '<xpress path>/include'[0m"
    echo "[31m[1mExample: ?> ./build-all.sh -DGLPK_PATH=/usr[0m"
    echo "[31m[1mPlease, install and choose at least one ILP Solver[0m"
    exit 1
else
    echo "Parameters: "$param
fi

cd grid
cmake CMakeLists.txt
cmake CMakeLists.txt
make
cd -

cd polygon
cmake CMakeLists.txt
cmake CMakeLists.txt
make
cd -

cd scp-solver
cmake CMakeLists.txt $param
cmake CMakeLists.txt $param
make
cd -

cd pre-solver
cmake CMakeLists.txt $param
cmake CMakeLists.txt $param
make
cd -

cd art-gallery-pg
cmake CMakeLists.txt $param
cmake CMakeLists.txt $param
make
cd -
