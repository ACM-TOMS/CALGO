#!/bin/bash

path=`dirname $0`
cd $path

echo ""
echo "*************************************************************"
echo "*           COMPILATION OF DELTAGAMMAINC MODULES            *"
echo "*************************************************************"
echo ""

gcc -w -O3 kernel.c Gfunc.c -lm -o Gfunc
if [ "$?" = "0" ]; then echo "  + compilation of module 'Gfunc': success"; else echo "  + compilation of module 'Gfunc': failure"; fi

gcc -w -O3 kernel.c deltagammainc.c -lm -o deltagammainc
if [ "$?" = "0" ]; then echo "  + compilation of module 'deltagammainc': success"; else echo "  + compilation of module 'deltagammainc': failure"; fi

gcc -w -O3 kernel.c check.c -lm -o check
if [ "$?" = "0" ]; then echo "  + compilation of module 'check': success"; else echo "  + compilation of module 'check': failure"; fi

echo ""
