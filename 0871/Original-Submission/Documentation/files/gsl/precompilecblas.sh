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
    elif [[ -z $cfile ]]
    then
      target=`basename $file ".h"`
      touch $targdir/log.txt
      #echo "Converting $dir/$file into $targdir/$file"
      echo "" >> $targdir/log.txt
      echo "$file" >> $targdir/log.txt
      /home/wschrep/precompiler/precompile $dir/$file $targdir/$file -x /home/wschrep/precompiler/mpieee-noconst.xml >> $targdir/log.txt 2>&1
      echo "/home/wschrep/precompiler/precompile $dir/$file $targdir/$file -x /home/wschrep/precompiler/mpieee-noconst.xml"
    elif [[ -z $header ]]
    then
      target=`basename $file ".c"`

      #echo "Converting $dir/$file into $targdir/$file"
      echo "" >> $targdir"/"log.txt
      echo "$file" >> $targdir"/"log.txt
      /home/wschrep/precompiler/precompile $dir/$file $targdir/$file -x /home/wschrep/precompiler/mpieee-noconst.xml >> $targdir/log.txt 2>&1
      echo "/home/wschrep/precompiler/precompile $dir/$file $targdir/$file -x /home/wschrep/precompiler/mpieee-noconst.xml"
    fi
  done
}

convertdir $1 $2

