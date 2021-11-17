File README.txt

This is a guide to archive XXX.zip attached to
Masao Kodama: "Algorithm XXX: A module for calculating cylindrical functions of 
complex order and complex argument," ACM Transactions on Mathematical Software.

This archive contains the following 6 files.

README.txt
File README.txt is this file and explains the contents of each file stored 
in this archive.

mod_zbes.f90
File mod_zbes.f95 stores MODULE mod_zbes describing the present algorithm.

examination.f90
This file stores PROGRAM examination, which is a sample program using MODULE 
mod_zbes. A method of implementation of this program is written in the 
beginning of this file.

examination_04.out
This is an outputted file from PROGRAM examination when kp=KIND(1.). 
kp is the kind type parameter in MODULE mod_zbes.

examination_08.out
This is an output file from PROGRAM examination when kp=KIND(1D0).

examination_16.out
This is an output file from PROGRAM examination when 
kp=SELECTED_REAL_KIND(20,550).

The above output files from PROGRAM examination.f90 is
obtained by using the following processor:
  OS: Windows Vista (64-bit version)
  CPU: Intel Core 2 Duo E4500
  Compiler: Intel Visual Fortran Compiler for Windows, Version 10.1.021

Author: Masao kodama
Address: 21-20, Gakuen 1 chome, Mastue-shi, Shimane-ken, 690-0825 Japan
Email: mkodama@mable.ne.jp

