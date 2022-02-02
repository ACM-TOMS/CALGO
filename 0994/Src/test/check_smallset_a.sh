# =============================================================================
# Check all matrices in "../mat_smallset_a" directory.
# =============================================================================
echo " "

DIR=../matrices/mat_smallset_a
BIN=../src/test_distance

# Get config.in
cp ../config.in .

# Remove previous results.
rm -f ${DIR}/*.out

# Process every matrix in "../mat_test" directory.
for i in $( ls ${DIR}/atest*.in )
do
  echo -n "Testing file: " $i
  fileOut=${i%.in}.out
  #### echo "Output file is: " $fileOut
  $BIN $i > $fileOut
  #### grep "Distance of input matrix:" $fileOut

  # Get the initial distance (the one in the name of the file).
  distInitial=$( echo $i | sed 's/^.*x//' | sed 's/\.in//' )
  #### echo "distInitial: " $distInitial

  # Get the computed distance by the code.
  distComputed=$( grep "Distance of input matrix:" $fileOut | \
                  sed 's/Distance of input matrix://' )
  #### echo "distComputed: " $distComputed

  # Check if everything is ok.
  if [ "$distInitial" -eq "$distComputed" ]
  then
    # Result was ok.
    resultOfTest="Ok"
    echo "    " $resultOfTest
  else
    # Failure.
    resultOfTest="!!! Failure !!!"
    echo "    " $resultOfTest
    #### echo "  "
  fi

  # Remove output file.
  #### rm $fileOut
done

# Remove config.in
rm config.in

echo

