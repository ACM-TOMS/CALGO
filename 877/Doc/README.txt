File README.txt

Author: Masao kodama
Address: 21-20, Gakuen 1 chome, Mastue-shi, 
         Shimane-ken, 690-0825 Japan
Email: mkodama@mable.ne.jp

This is a guide to the algorithm attached to
Masao Kodama: "Algorithm XXX: A subroutine package for Bessel functions of 
complex order and a nonnegative argument," ACM Transactions on Mathematical 
Software.

Archive file soft_bes.gz contains the following nine files.
Integer kp below denotes the kind type parameter used in the present algorithm.

README.txt
File README.txt is this file, and explains the contents of each file
included in archive file soft_bes.gz.

pack_bes.f95
File pack_bes.f95 contains the present algorithm.

test_program.f95
File test_program.f95 includes PROGRAM test_bes, which checks the present 
algorithm. PROGRAM test_bes contains Tests 1-9.

test_program_04.out
One of the output files from PROGRAM test_program. 
test_bes_04.out is the output file when kp=KIND(1.0).
The processor used here is as follows.
OS: Windows XP
CPU: Athlon(tm) XP 2200+
Compiler: Fortran&C Package

test_program_08.out
One of the output files from PROGRAM test_program. 
The file test_bes_08.out is the output file when kp=KIND(1D0).
The processor used here is as follows.
OS: Windows XP
CPU: Athlon(tm) XP 2200+
Compiler: Fortran&C Package

test_program_16.out
One of the output files from PROGRAM test_program. 
The file test_bes_16.out is the output file when kp=SELECTED_REAL_KIND(20,550).
The processor used here is as follows.
OS: Windows XP
CPU: Athlon(tm) XP 2200+
Compiler: Fortran&C Package

test_coulcc.f95
The file test_coulcc.f95 has PROGRAM test_coulcc, which is a test program
for Algorithm COULCC.
The test program contains Test 10.

test_coulcc.out
The file test_coulcc.out is the output file from PROGRAM test_coulcc when 
kp=KIND(1D0).
The processor used here is as follows.
OS: Windows XP
CPU: Athlon(tm) XP 2200+
Compiler: Fortran&C Package

transformer.f95
File transformer.f95 includes PROGRAM transformer, which transforms programs
from Fortran 77 to Fortran 95.

