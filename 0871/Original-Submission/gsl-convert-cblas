#!/bin/bash


convertdir(){
  local dir=$1
  local targdir=$2
  echo "dir=$dir, target dir=$targdir"
  
  
  for file in `ls -1 $dir`;
  do
    check=$dir"/"$file
    if [[ -d $check ]]
    then
      local newdir=$dir"/"$file
      local newtarg=$targdir"/"$file
      convertdir $newdir $newtarg
#    else
#      echo "$check is not a dir"
    fi
  
    local target=""
    local header=`echo $file | grep "\.h"`
    local cfile=`echo $file | grep "\.c"`
    local both=`echo "$header$cfile" | grep -v "_pre"`
  
    if [[ -z $both ]]
    then
      #skip file
      echo -n 
    else
      target=`basename $file ".h"`
      touch $targdir/log.txt
      #echo "Converting $dir/$file into $targdir/$file"
      echo "" >> $targdir/log.txt
      echo "$file" >> $targdir/log.txt
      cat $dir/$file | sed s/"#define BASE double"/"#define BASE MpIeee"/ | \
       sed s/"#define BASE float"/"#define BASE MpIeee"/ | \
       sed s/"#define ACC_TYPE  float"/"#define ACC_TYPE MpIeee"/ | \
       sed s/"#define ACC_TYPE  double"/"#define ACC_TYPE MpIeee"/  > $targdir/$file
        
      echo "converting BASE double to MpIeee in ".$file
      
    fi
  done
}

if(( $#!=2 ));
then
	echo "USAGE: $0 <gsl cblas original> <gsl cblas precompiled>";
	exit 1;
fi

convertdir $1 $2

