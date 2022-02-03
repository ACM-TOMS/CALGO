#!/bin/bash

counter=1;

function testcount()
{
  echo -n "$counter."
  if(( $counter < 10 )); 
  then
    echo -n " "; 
  fi;
  (( counter=$counter+1 ))
}

#########################
# function test
#
# run argument given and put output
# in a .new file, then compare the .new and .out
# if diff gives back results -> test failed!
function test()
{
  testcount
    
  echo -n " testing : $1.cpp  ->  ";
  ../precompile $1.cpp $1_pre.cpp -x convert.xml 2>/dev/null >/dev/null;

  local ok;
  ok=`diff $1_pre.cpp $1_good.cpp`;
  
  if [[ $ok  == "" ]] ;
  then
    echo "ok";
    rm -f $1_pre.cpp;
  else
    echo "TEST FAILED!!!";
  fi;
}


#########################
# main
#

rm -f *_pre* *~

test "test"       
test "test2"      
test "test3"      
test "test4"      
test "array"      
test "functest"   
test "funcdeclare"


echo -n "testing : functest.cpp with skip.conf -> ";
../precompile functest.cpp functest_pre2.cpp -x convert.xml -c skip.conf  2>/dev/null >/dev/null;
ok=`diff functest_pre2.cpp functest_good2.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "ok";
else
  echo "TEST FAILED!";
fi;



echo -n "testing : functest.cpp with skip2.conf -> ";
../precompile functest.cpp functest_pre3.cpp -x convert.xml -c skip2.conf 2>/dev/null >/dev/null;
ok=`diff functest_pre3.cpp functest_good3.cpp`;
if [[ $ok  == "" ]] ;
then
  echo "ok";
else
  echo "TEST FAILED!";
fi;


