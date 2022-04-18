
An exact algorithm for the three-dimensional bin packing problem
================================================================

The code consists of two parts:

  3dbpp.c    : The callable C-code which solves a three-dimensional
               bin packing problem.

  test3dbpp.c: A main algorithm which either reads input from a file
               or generates 10 instances of a given size/type, and 
               calls 3dbpp.c for the solution. 

To compile the code use one of the following commands

  gcc -ansi -o 3dbpp -O5 3dbpp.c test3dbpp.c -lm  (gnu C)
  cc  -Aa   -o 3dbpp +O4 3dbpp.c test3dbpp.c -lm  (HP-UX C)

The generated executable file "3dbpp" can either read an instance from 
a file or randomly generate 10 instances.

If the instance is read from a file, five arguments should be given:
  filename  A filename in which the test instance is found. The format 
            of the file is:
               n W H D
               w_1 h_1 d_1
               :
               w_n h_n d_n
            where 
               n is the number of items, 
               W,H,D is the size of the bin, and 
               w_j,h_j,d_j is the size of box j.
  nodelimit maximum number of decision nodes to be explored in the
            main branching tree. If set to zero, the algorithm will
            run until an optimal solution is found (or timelimit or
            iterlimit is reached). Measured in thousands (see IUNIT).
  iterlimit maximum number of iterations in the ONEBIN algorithm
            which packs a single bin. If set to zero, the algorithm will
            run until an optimal solution is found (or timelimit or
            nodelimit is reached). Measured in thousands (see IUNIT).
  timelimit Time limit for solving the problem expressed in seconds.
            If set to zero, the algorithm will run until an optimal
            solution is found; otherwise it terminates after timelimit
            seconds with a heuristic solution. 
  nodeused  returns the number of branch-and-bound nodes investigated,
            measured in thousands (see constant IUNIT in 3dbpp.c).
  iterused  returns the number of iterations in ONEBIN algorithm,
            measured in thousands (see constant IUNIT in 3dbpp.c).
  timeused  returns the time used in miliseconds
  packingtype 
            Desired packing type. If set to zero, the algorithm will
            search for an optimal general packing; if set to one, it
            will search for a robot packing.
            will search for a robot packing.
If the code should randomly generate 10 instances, seven arguments 
should be given:
  n         The size of the instance, i.e., number of boxes.
  bindim    The size of the bin, typically 40-100.
  type      An integer saying which randomly generated instance should
            be generated. A value between 1-9 selects one of the instance 
            types described in the above papers. Value 10-11 generates
            1D and 2D instances.
  nodelimit as above
  iterlimit as above
  timelimit as above
  packingtype
            as above

Results are written to standard output.

Thus for instance choosing n=30, bindim=100, type=1, nodelimit=0,
iterlimit=0, timelimit=0, packingtype=1 one should get
the following output:

3DBPP PROBLEM 30 100 1 0 0 0 1
 1 : lb  9 z  9 node         0 iter         3 time   0.00
 2 : lb 11 z 11 node        41 iter       438 time   0.32
 3 : lb  8 z  8 node       239 iter     12392 time  12.79
 4 : lb  7 z  7 node         0 iter       978 time   1.14
 5 : lb 10 z 10 node      2194 iter     91653 time  79.18
 6 : lb  7 z  7 node         0 iter         2 time   0.00
 7 : lb  8 z  8 node      1716 iter     30606 time  32.02
 8 : lb  9 z  9 node         0 iter         6 time   0.00
 9 : lb  9 z  9 node         0 iter        38 time   0.02
10 : lb  8 z  8 node         5 iter       393 time   0.27
as well as some average values on the above values.

"lb" is the lower bound, "z" is the objective of the found solution,
"node" is the number of branch-and-bound nodes in thousands, and
"iter" is the number of iterations in the ONEBIN algorithm in thousands,
and "time" is the CPU-time used for solving the problem.

(c) Copyright 1998, 2003, 2005

  David Pisinger                        Silvano Martello, Daniele Vigo
  DIKU, University of Copenhagen        DEIS, University of Bologna
  Universitetsparken 1                  Viale Risorgimento 2
  Copenhagen, Denmark                   Bologna, Italy

This code can be used free of charge for research and academic purposes 
only. 

