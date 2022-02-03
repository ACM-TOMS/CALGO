#!/bin/bash

ok=`rm -f *_pre* *~`

echo -n "testing : test.cpp  -> ";
../precompile test.cpp test_pre.cpp -x convert.xml ;
ok=`diff test_pre.cpp test_good.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;


echo -n "testing : test2.cpp -> ";
../precompile test2.cpp test2_pre.cpp -x convert.xml  ;
ok=`diff test2_pre.cpp test2_good.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;


echo -n "testing : test3.cpp -> ";
../precompile test3.cpp test3_pre.cpp -x convert.xml ;
ok=`diff test3_pre.cpp test3_good.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;


echo -n "testing : test4.cpp -> ";
../precompile test4.cpp test4_pre.cpp -x convert.xml ;
ok=`diff test4_pre.cpp test4_good.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;


echo -n "testing : array.cpp -> ";
../precompile array.cpp array_pre.cpp -x convert.xml ;
ok=`diff array_pre.cpp array_good.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;


echo -n "testing : functest.cpp -> ";
../precompile functest.cpp functest_pre.cpp -x convert.xml ;
ok=`diff functest_pre.cpp functest_good.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;



echo -n "testing : functest.cpp with skip.conf -> ";
../precompile functest.cpp functest_pre2.cpp -x convert.xml -c skip.conf
ok=`diff functest_pre2.cpp functest_good2.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;



echo -n "testing : functest.cpp with skip2.conf -> ";
../precompile functest.cpp functest_pre3.cpp -x convert.xml -c skip2.conf
ok=`diff functest_pre3.cpp functest_good3.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;



echo -n "testing : funcdeclare.cpp -> ";
../precompile funcdeclare.cpp funcdeclare_pre.cpp -x convert.xml
ok=`diff funcdeclare_pre.cpp funcdeclare_good.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;

echo -n "testing : test_operation_extension.cpp -> ";
../precompile test_operation_extension.cpp test_operation_extension_pre.cpp -x operationtest.xml
ok=`diff test_operation_extension_pre.cpp test_operation_extension_good.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "TEST OK";
else
  echo "TEST FAILED!";
fi;


