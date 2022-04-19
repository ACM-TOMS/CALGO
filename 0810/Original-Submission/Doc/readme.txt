
 
                   SLEIGN2: THE README.TXT FILE
  
   01 March 2001; P.B. BAILEY, W.N. EVERITT AND A. ZETTL
  
  PLEASE read carefully this readme.txt file and the intro.tex file before 
 using the SLEIGN2 package for the first time.
  
 	There are eleven files in the SLEIGN2 package as follows:
 
 	Two ASCII files:
 	autoinput.txt  readme.txt
 
 	Six FORTRAN files:
 	coupdr.f  drive.f  makepqw.f  sepdr.f  sleign2.f  xamples.f
 
 	Three AMS-LaTeX files (these files can be compiled in the UNIX latex
 compiler, and then printed out in hard copy):
 	help.tex  intro.tex  xamples.tex   
  
   To run one of the examples in the xamples.f file in an UNIX environment 
 with a Fortran77 (or Fortran77) compiler, enter the following command:
  
 f77 xamples.f  drive.f  sleign2.f  -o  xamples.x
 (Replace f77 by f90 if you want to use the Fortran90 compiler.)
  
   This will create the executable file xamples.x and object files drive.o
 and sleign2.o. Now run xamples.x whenever you want to work an example from 
 the list of examples in the xamples.f file. The hard printed copy of the
 file xamples.tex provides detailed information on each one of the 32
 chosen examples in the file xamples.x.
  
   To run your own problem proceed as follows:
  Step 1: Enter the command 
        f77  makepqw.f  -o  makepqw.x
  Step 2: Run
        makepqw.x  
   (This interactive program will ask you for a file name - this must end 
 in .f - for example, if your chosen problem has the name bloggs then enter 
 bloggs.f). This file, after makepqw.x has been run, will contain the 
 subroutines for p,q,w and, if necessary, the functions u, v and U, V which 
 are used to define singular limit-circle boundary conditions at one or both 
 endpoints.
  Step 3: Enter the command
        f77  bloggs.f  drive.f  sleign2.f  -o  bloggs.x
   (drive.f and sleign2.f can be replaced with drive.o and sleign2.o 
 if these .o files are available to speed up the compilation; these -o 
 files are created by the first compilation.)
  Step 4: Run
        bloggs.x  
   (The user is asked to provide information for the code to run
 bloggs.x, for example: boundary conditions, eigenvalue indexes, 
 numerical tolerances, name for report file if desired. See the autoinput.txt 
 file for instructions on how to automate the input of the required 
 information).
  
  The file sepdr.f is a sample driver for separated regular and singular
 boundary conditions; coupdr.f is a sample driver program for coupled regular
 and singular boundary conditions. Experienced users who want to bypass the 
 extensive user friendly interface provided in drive.f and makepqw.f, and
 use their own driver may wish to look at these two sample drivers.

  Note that the above procedures may have to be modified for non UNIX 
 environments, e.g. DOS or APPLE.
 
  Note also that in running xamples.x or bloggs.x the user has access to
 the interactive help device; at any point where the program halts  
 type h <ENTER> for access; to return to the point in the program where
 help was accessed, type r <ENTER>. Additional information on the help 
 device can be found in the intro.tex file. The whole of the help data
 can be printed out separately from the file help.tex and the user is
 advised to have a copy available for consultation when working with the
 code files for the first time.
 
  See the autoinput.txt file for instructions on how to bypass these 
 program halts in xamples.x and bloggs.x
  
  All six of the FORTRAN .f files are supplied in single precision; 
 to convert these files to double precision replace the string `  REAL' by 
 the string `  DOUBLE PRECISION' throughout. (Note the two spaces in front of 
 REAL and in front of DOUBLE PRECISION - these spaces are important.)  
 In UNIX this can be effected within the vi editor as follows:
  
 :1,$ s/  REAL/  DOUBLE PRECISION
  
  It is recommended that the user try the program in single precision and 
 switch to double precision as required.
  
  Additional information on the SLEIGN2 package can be found in the intro.tex
 and help.tex files. See also a hard printed copy of the file bailey.tex,
 available in the Zettl/Sleign2 directory detailed below, which contains
 an account of the analytical and numerical properties of the SLEIGN2
 code, together with an extensive list of references.
  
  All of the eleven files in the SLEIGN2 package can be obtained by anonymous
 ftp from the computer
        ftp.math.niu.edu
 and the directory
        /pub/papers/Zettl/Sleign2
  
  A sample session using the traditional UNIX or DOS ftp client is as 
 follows:
  
  Step 1. At the UNIX or DOS prompt type
              fop  ftp.math.niu.edu
  
  Step 2. In response to the request for user id type
              fop
  
  Step 3. In response to the request for username type your e-mail address;
             for example:
             w.n.everitt@bham.ac.uk
  
  Step 4. At the ftp> prompt enter:
             cd  pub/papers/Zettl/Sleign2
             get  readme.txt
             quit
  
  This will transfer the readme.txt file to your current directory and close 
 the connection.
  
  At the ftp> prompt you can also enter 'ls' to see the listing of other 
 files available in the Sleign2 directory.
  
  Users who prefer World Wide Web software such as Netscape 
 or lynx can specify the URL
            ftp://ftp.math.niu.edu/pub/papers/Zettl/Sleign2
 to access the Sleign2 directory; this software will then show the list 
 of files available in it.
 
  Replacing the directory "Sleign2" with the directory "Pub papers" and
 repeating the above procedures will make available a number of recent papers
 which are related to the SLEIGN2 package. For example, the paper "BEWZ" 
 (Bailey, Everitt, Weidmann and Zettl) contains some results which can be 
 combined with SLEIGN2 to approximate the continuous (essiential) spectrum 
 of singular limit-point Sturm-Liouville problems.
  
  All eleven files of the SLEIGN2 package and a number of recent publications
 related to it, can also be accessed through the web page:
           http://www.math.niu.edu/~zettl/SL2/
  
  All suggestions, comments and criticisms are welcome; please send all 
 comments to Tony Zettl at
          zettl@math.niu.edu
  
     Paul Bailey, Norrie Everitt and Tony Zettl 
 (with the assistance of Burt Garbow.)
  
  Acknowledgement. The authors are grateful to their colleagues Eric Behr, 
 Qingkai Kong and Hongyou Wu for help and advice at a number of stages 
 in the development of this program. Some of the theoretical underpinnings 
 of the algorithm for coupled boundary conditions were obtained jointly 
 with Michael Eastham, Qingkai Kong and Hongyou Wu.

  A special thanks to Eric Behr for his help throughout the development of
 the code, for setting up the public access through the Internet and the
 world wide web, and for informed advice.

 Department of Mathematical Sciences, Northern Illinois University 
 DeKalb, IL 60115-2888, USA



