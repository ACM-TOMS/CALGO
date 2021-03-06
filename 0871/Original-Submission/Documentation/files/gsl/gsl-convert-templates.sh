#!/bin/bash


convertFunc(){
  local origfile=$1
  local targfile=$2
  
  cat $origfile | sed s/"#define BASE double"/"#define BASE MpIeee"/        | \
       sed s/"#define BASE long double"/"#define BASE MpIeee"/   | \
       sed s/"#define ATOMIC double"/"#define BASE MpIeee"/      | \
       sed s/"#define ATOMIC long double"/"#define BASE MpIeee"/ | \
       sed s/"#define BASE float"/"#define BASE MpIeee"/         | \
       sed s/"#define ATOMIC float"/"#define BASE MpIeee"/ > $targfile
       
  echo "converting BASE and ATOMIC to MpIeee in $targfile"
      
}

if(( $#!=2 ));
then
	echo "USAGE: $0 <original template file> <converted template file>";
	exit 1;
fi

convertFunc $1 $2


