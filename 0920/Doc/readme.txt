
SFSDP
(A Sparse Version of Full SemiDefinite Programming Relaxation for 
Sensor Network Localization Problems)

Sunyoung Kim, Masakazu Kojima, Hayato Waki and Makoto Yamashita

----------------------------------------
Index

1. Overview
2. System Requirements
3. Installation Guide
4. Usage of SFSDP
5. Acknowledgements
6. License
7. E-mail address
8. History

----------------------------------------

1. Overview

SFSDP is a Matlab package for solving sensor network localization 
problems. The package contains  four functions, SFSDP.m, SFSDPplus.m, 
generateProblem.m, and test_SFDP.m. Numerical examples are included. 

The function SFSDP.m is a Matlab implementation of the semidefinite
programming relaxation proposed in the recent paper by Kim, Kojima and
Waki for sensor network localization problems, as a sparse version of
the full semidefinite programming relaxation (FSDP) by Biswas and
Ye. To improve the efficiency of FSDP, SFSDP.m exploits the aggregated
and correlative sparsity of a sensor network localization problem.


2. System Requirements

The following softwares are required for SFSDP.

  i. MATLAB R2006b or later.
	-- available from Mathworks inc.

  ii. sdpa.7.3.1 or later  
	-- available from http://sdpa.indsys.chuo-u.ac.jp/sdpa/
	to call sedumiwrap.m from SFSDP for solving an SDP relaxation problem
  and/or 
  iii. SeDuMi 1.1R3 or later
	-- available from http://sedumi.ie.lehigh.edu/ 
	to call SeDuMi from SFSDP for solving an SDP relaxation problem

Note: (1) We recommend to use sdpa for solving an SDP relaxation  
because sdpa is faster than SeDuMi for large-scale sensor network 
localization problems. See Section 5 in User Manual of SFSDP 
for the numerical results. 

(2) In the case of ii, add a path to sedumiwrap,  
       ~/sdpa/sdpa.7.3.1/mex/
    in Matlab.

3. Installation Guide

  1 Download SFSDP from:
	http://www.is.titech.ac.jp/~kojima/SFSDP

  2 Unpack it by either:
	double-clicking the icon, 
	or typing on the Terminal: $ tar zxvf SFSDP.tar.gz

	You will get a folder named SFSDP.

  3 Move the folder anywhere you like, for example, to
	/Users/smith/matlab/SFSDP (if you are smith),
	
  4 Invoke Matlab, and cd to the folder of SFSDP.
	>> cd /Users/smith/matlab/SFSDP

  5 Set a path to SFSDP. For example,
	Choose: File->Set Path
	Push 'Add with Subfolders...' button, and 
	Select the folder SFSDP.
	Finally, Save it by pushing 'Save' button.

  6 Now you are ready to use SFSDP.

NOTE: If you want to use Matlab in the mode 'matlab -nodisply', 
it is useful to add the path of the directory SFSDP in 
your startup.m. If your SFSDP is installed in 
/Users/smith/matlab/SFSDP, add the following command in your 
startup.m:

addpath(genpath('/Users/smith/matlab/SFSDP'));


4. Usage of SFSDP

We assume that SDPA and/or SeDuMi have been already installed in 
your computer and that the paths of their folders and subfolders 
have been already added in your Matlab search path.

The usage of SFSDP is explained in the user manual available at
http://www.is.titech.ac.jp/~kojima/SFSDP


5. Acknowledgments

The authors are grateful to Professor Yinyu Ye for the original
version of FSDP, and Professor Kim-Chuan Toh for Matlab programs 
refineposition.m and procrustes.m, and valuable suggestions.
The Matlab programs of Professor Kim-Chuan Toh are distributed
under the GNU General Public License 2.0 and available at
http://www.math.nus.edu.sg/~mattohkc/SNLSDP.html


6. License

This software is distributed under
GNU General Public License 2.0

7. E-mail address

kojima.m.aa-sfsdp@m.titech.ac.jp

Please send a message if you have any questions, ideas, or 
bug reports.


8. History

Version 1.11 (July 31, 2009)
(1) Users can use SDPA instead of SeDuMi for solving an SDP 
relaxation problem.  See `readmeV111.txt' for more detail.

Version 1.01 (August 03, 2008)
(1) A release after tuning the first version.
(2) Became available at http:www.is.titech.ac.jp/~kojima/SFSDP.

Version 1.00 (July 28, 2008)
(1) The first version. 

 
