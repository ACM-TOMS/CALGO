program checktmessy
! On Windows, run the batch file checktmessy.bat that compares
! test output with the distribution files.
! The results of the comparison are in d_summary.
   call system("echo off") ! avoid printing unwanted info
! Fetch result from tmessy; no prompting for overwrite.
   call system("copy ..\tmessy\d_newresult . /y > nul")
! Compare results from distribution file found in \result_orig.   
   call system("checktmessy > d_summary")
! Show summary of comparison.
   call system("type d_summary")
end program