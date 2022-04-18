
Instructions


The library routines are in the file lib.f. These routines take a
compressed reference and place its components in the pair of COMMON
blocks CREFNO and CREFST.

There are two example programs 

Example 1 -- Building a bibtex database
---------------------------------------

Edit the routines SETUP in shared.f and OUTSET in bibeg.f
to ensure the OPEN statements are valid.
 
Compile and link lib.f shared.f bibeg.f using your favourite
Fortran 77 compiler.

Running the resultant object code should produce the bibtex
database in the file specified in OUTSET.

Example 2 -- Cumulative index sorted on SHARE index
---------------------------------------------------

Edit the routines SETUP in shared.f and OUTCUM in cumeg.f
to ensure the OPEN statements are valid.
 
Compile and link lib.f shared.f cumeg.f using your favourite
Fortran 77 compiler.
NOTE: If you are running on a UNIX system you may need to use
the file cumeg.unix.f because of the use of backslashes in
write statements.

Running the resultant object code should produce the LaTeX
input necessary to generate a cumulative index. 

