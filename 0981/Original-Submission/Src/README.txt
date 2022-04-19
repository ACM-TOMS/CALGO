This folder contains:

    ./TalbotSuiteDE: subfolder for the Talbot Suite DE software.
    ./ex*          : subfolders with sequential (./1SEQ) and
                     parallel (./2PAR) examples.

See the User Guide (TsuiteDE_UserGuide.pdf) for details.
See directoryTree.pdf for the folder tree.


To build and execute an executable related to a particular example,
we provide:
  - runme.sh for Linux (remember to make executable the file)
  - runme.bat for Windows
  - runme.m for MATLAB


In addition, every example contains a file ("output.txt") with the
corresponding output results. When there are two main programs, for
example "main_ACCURACY" and "main_TIMES", there are two text files
containing their results: "output_ACCURACY.txt" and "output_TIMES.txt"
respectively.


See the User Guide for important remarks about the solution to some
possible errors issued during the compilation or execution of mixed
C-OpenMP/MATLAB code under Windows and under Linux.

