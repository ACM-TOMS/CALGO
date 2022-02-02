AGSol (Art Gallery Solver)      v1.0.2      16/Dec/2015
===============================================================================
  Code Authors:
    Davi Colli Tozoni - davi.tozoni@gmail.com
    Marcelo Castilho Couto - coutomarcelo@gmail.com
  
  Concept and Design: 
    Davi Colli Tozoni - davi.tozoni@gmail.com
    Marcelo Castilho Couto - coutomarcelo@gmail.com
    Pedro Jussieu de Rezende - rezende@ic.unicamp.br
    Cid Carvalho de Souza - cid@ic.unicamp.br

This package contains a software capable of optimally solving the Art Gallery
Problem (AGP), one interesting NP-hard problem from the Computational Geometry
field. The algorithm implemented in this solution, which can be today
considered the state-of-the-art technique on the AGP, can be found in details
in the following paper:

[1] Davi C. Tozoni, Pedro J. de Rezende, Cid C. de Souza, 2013. A Practical
Iterative Algorithm for the Art Gallery Problem using Integer Linear
Programming. 

All experiments displayed in [1] were performed employing AGSol. 

In this README file, our objective is to present general information about our
solution. First, we introduce the Art Gallery Problem and talk about some of
its fundamental characteristics. After that, the algorithm in [1] is briefly
presented and discussed. In section 3, we focus on our implementation and
explain the main modules of AGSol. The subsequent sections describe the
necessary procedure for installing AGSol and how to run it in your machine.
This description also includes the possible inputs and outputs of AGSol.
Finally, our README ends by showing simple examples on how to properly use
AGSol.

***IMPORTANT: This software was initially designed to be use in a Linux OS.
Attempts to install and run it on a different platform may require a higher
level of proficiency on the part of the user.

-------------------------------------------------------------------------------


1) The Art Gallery Problem (AGP) 

The problem tackled by this code, known as the Art Gallery Problem (AGP), is a
NP-Hard problem that consists in finding the minimum number of guards that are
sufficient to ensure that an art gallery (represented by a simple polygon) is
fully guarded, assuming that a guard’s field of view covers 360 degrees as well
as an unbounded distance. Observe the example below, where 2 guards (G1 and G2)
are sufficient to cover the whole polygon:
 ________________
|                \
|   G1            \
|                  \ 
|             _____|
|            |
|            |________
|                     |
|       G2            |
|_____________________|


The polygons to be covered can also have holes. In an art gallery, a hole can 
be seen as a pillar or any other obstacle to the guard's vision. See the 
following example:
   __________________
  /                 /____________
 /      G1                       \
|                      __         |
|           /\        |  |        |
|          /  \       |__|        |
|         |    \                 /
|         |_____|               /
|_                        G2   |
  |____________________________|
  
  
To present a motivation for solving the problem, imagine that the National
Museum of Rio de Janeiro, in Brazil, is promoting an exhibition in which new
archaeological discoveries will be presented. Imagine also that the owners of
the exhibition want to install cameras in the building in such a way the whole
interior is guarded, using, however, a limited budget for this task. Thus, how
is the best way to position the cameras? An answer for this problem could be
found with AGSol by reducing the input gallery into a polygon.

The natural motivation just described represents the original Art Gallery
Problem, in which guards can be positioned anywhere in the polygon and where
the objective is to watch over the interior of the gallery. There are, however,
several other variants of the problem, such as the Art Gallery Problem With
Fixed Guard Candidates (AGPFC), where a viable solution consists of a set
guards that guarantee surveillance of the entire polygon, while having their
placement restricted to a finite number of positions. In contrast, consider the
Art Gallery Problem With Witnesses (AGPW) that asks for a set of guards able to
observe merely a given finite set of points inside the polygon, while guard
placement may be unrestricted. 


1.a) Algorithmic achievements

During the last years, there were several attempts to solve the Art Gallery
Problem. Because the AGP is a NP-hard problem, many of these attempts involved
the use of heuristics and approximation algorithms. Quite recently, however,
some authors have focused on the search for proven optimal solutions to the
problem.  The timeline below presents some of the major recent works on the
AGP. See that [1], as previously said, corresponds to the algorithm implemented
in AGSol.
 
[2]          [3,4]      [5]       [6]       [7]       [8]                 [1]
 |-----...-----|---------|---------|---------|---------|---------|---------|
1987          2007      2008      2009      2010      2011      2012      2013


Works:

[2] S. K. Ghosh. Approximation algorithms for art gallery problems. In Proc.
Canadian Inform. Process. Soc. Congress, pages 429–434, Mississauga, Ontario,
Canada, 1987.  Canadian Information Processing Society.

[3] M. C. Couto, C. C. de Souza, and P. J. de Rezende. An exact and efficient
algorithm for the orthogonal art gallery problem. In Proc. of the XX Brazilian
Symp. on Comp.  Graphics and Image Processing, pages 87–94. IEEE Computer
Society, 2007.

[4] Y. Amit, J. S. B. Mitchell, and E. Packer. Locating guards for visibility
coverage of polygons. In ALENEX, pages 1–15, New Orleans, Lousiana, January
2007. SIAM.

[5] A. Bottino and A. Laurentini. A nearly optimal sensor placement algorithm
for boundary coverage. Pattern Recognition, 41(11):3343–3355, 2008.

[6] M. C. Couto, P. J. de Rezende, and C. C. de Souza. An exact algorithm for
an art gallery problem. Technical Report IC-09-46, Institute of Computing,
University of Campinas, Nov. 2009.

[7] T. Baumgartner, S. P. Fekete, A. Kröller, and C. Schmidt. Exact solutions
and bounds for general art gallery problems. In Proceedings of the SIAM-ACM
Workshop on Algorithm Engineering and Experiments, ALENEX 2010, pages 11–22.
SIAM, 2010.

[8] A. Bottino and A. Laurentini. A nearly optimal algorithm for covering the
interior of an art gallery. Pattern Recognition, 44(5):1048–1056, 2011.


1.b) Experiments:

In our experiments, whose results are fully reported in [1], AGSol was tested
and compared against works [7] and [8] and achieved a significant improvement
over the previously published techniques. Using AGSol, optimal solutions are
reached in a few minutes for almost the totality of the publicly available
instances. Not only AGSol was able to solve instances with up to 2500 vertices,
something without precedent in the literature, but it also proved to be quite
robust, achieving optimality in more than 98% of the 2400+ instances tested.

The classes of polygons considered in our experiments are the following:

+ Simple: this class refers to random simple polygons (without holes).
+ Orthogonal: this class comprises random polygons generated respecting the
property that all edges are parallel to the x-axis or y-axis.
+ von Koch: class of polygons created in Couto et al.’s work [6] based on a
modified version of the von Koch curve.
+ Simple-simple: class that comprises the instances where the boundary and all
the holes are simple polygons.
+ Ortho-ortho: class that comprises the instances where the boundary and all
the holes are orthogonal polygons.

Obs.: all polygons experimented in [1] are located in the 'instances' folder of
this package. 

See below a summary of the run time results of our software when solving
different classes of polygons with 100 vertices. 

   Instance class   |     Time (sec)
++++++++++++++++++++++++++++++++++++++++++
     Simple         |       35.04
------------------------------------------                    
    Orthogonal      |       20.64
------------------------------------------
     von Koch       |       76.45
------------------------------------------
  Simple w/ holes   |      498.69
------------------------------------------
Orthogonal w/ holes |      433.03

-------------------------------------------------------------------------------


2) The Software Package

This program is divided into 5 main modules. As presented in the package, each
of these modules is in a separate folder:


2.a) art-gallery-pg: 
    Represents the main algorithm. See the sketch below: 

AGP Algorithm
------------------------------------------------------------------------------
1:   Set UB ← |V| and LB ← 0
2:   Select the initial witness set D
3:   Select the initial guard candidate set C ⊇ V
4:   while (UB = LB) and (MAXTIME not reached) do
5:      Solve AGPW(D), set Gw ← optimal solution of AGPW(D) and 
6:       LB ← max{LB, |Gw|}
7:      Solve AGPFC(C), set Gf ← optimal solution of AGPFC(C and 
8:       UB ← min{UB, |Gf|}
9:      if (UB = LB) then
10:         Update D and C
11:     end if
12:  end while
13:  return Gf

For more information, see Section 4 of [1].


2.b) grid: 
    Responsible for constructing the arrangement and identifying Light AVPs.

    For the construction of the arrangement, we implemented a class that
extends the class CGAL::Arr_observer<Arrangement> of CGAL. This class defines
some actions that can be taken when events are thrown in the arrangement
construction, like the creation of a new edge or an existing edge splitting
into two.

    In its implementation, CGAL uses the famous DCEL structure to save
information on the current arrangement. For more details on how the
implementation of arrangements is done in CGAL, see the following web page:

[9] CGAL. 2015. Computational Geometry Algorithms Library. (2015). www.cgal.org 
(last access January 2015).

    The complexity of the resulting arrangement is of great importance for the
solver’s performance. For a generic set S of points in a polygon P, we have
O(|S| · |P|) visibility edges. Because of this, the complexity of constructing
the visibility arrangement of S, as well as the number of AVPs in the
arrangement, has in the worst case, considering that all of them intersect, a
final complexity of O(|S|^2 · |P|^2).
    

2.c) polygon: 
    Responsible for some important geometric operations, like the visibility
computation. Here, the Polygon_2 and Polygon_with_holes_2 classes from CGAL
were extended to include new and custom functions.
    
    The visibility computation algorithm implemented can be found in the
following paper: 

[10] B. Joe and R. B. Simpson. Visibility of a simple polygon from a point.
Report CS-85-38, Dept. Math. Comput. Sci., Drexel Univ., Philadelphia, PA,
1985.


2.d) pre-solver: 
    Corresponds to our "black-box" for solving ILP (SCP) instances. Here you
can find different modes for solving them. The current options are:
    1- GLPK + Lagrangian Heuristic
    2- GLPK
    3- XPRESS + Lagrangian Heuristic
    4- XPRESS

    In this module, we also implemented a method for removing redundant rows
and columns from the original SCP, in order to improve the performance of our
solution.


2.e) scp-solver:
    Contains several libraries, one for each different way of solving an SCP
instance. The "black-box" implemented in 'pre-solver' folder calls methods from
one or more of the libraries implemented here. 
    
    As an example, If the software is currently solving the AGP using GLPK,
then the pre-solver will be making calls to the SolverPLIGlpk.so library. 
    
    One of the possible methods to be used is the Lagrangian Heuristic. 
The heuristic implemented is based on the following work by Beasley: 

[11] John E. Beasley. 1993. Lagrangian Relaxation. In Modern Heuristic
Techniques for Combinatorial Problems, Colin R. Reeves (Ed.). John Wiley &
Sons, Inc., New York, NY, USA, 243–303.


Obs.: If you want to use this software with another ILP Solver package, it will
be necessary to write your own SolverPLIXXX.[C|h]. Moreover, it will also be
necessary to include a call to your new library in 'pre-solver' package
(PreSolver.C).

-------------------------------------------------------------------------------

3) How to install

    First of all, it is very important to notice that this software was
designed to be executed in a GNU/Linux system. Attempts to install and run it
on a different platform may require a higher level of proficiency on the part
of the user.

    In order to install it correctly, it is necessary to have an installed
version of CGAL in your machine, preferably the same used in our experiments
(version 3.9). This CGAL version can be downloaded with the following link:
    - CGAL version 3.9: http://gforge.inria.fr/frs/download.php/file/29125/CGAL-3.9.tar.gz

    In addition, you will also need to have at least one ILP solver
installed in your machine. This solver can be:
    - GLPK version 4.52: http://ftp.gnu.org/gnu/glpk/glpk-4.52.1.tar.gz
    - XPRESS version 7.0

    It does not mean that employing a different version of these libraries will
cause the program to fail. However, it is safer to use the mentioned versions,
since all experiments were done using them.

    Below you find a list of other software solutions known to be necessary in
order to correctly install AGSol:
    - cmake (2.8.7)
    - libgmp (libgmp.so.10.2)
    - libmpfr (libmpfr.4.1.0)
    - g++ (4.6.3)
    - libboost (including libboost-thread) (libboost 1.46.1)

    It is important to notice that these auxiliary libraries must be installed
before CGAL installation.

    Finally, to install the AGP Solver, you must run the command below in the
root directory of the package:

    ?> ./build-all.sh -DGLPK_PATH=<path in the system> -DXPRESS_PATH=<path in
     the system>

    Where the <path in the system> expression means that the solver's library
is located in '<path in the system>/lib' and that the header folder of the
solver can be found in '<path in the system>/include'.

    If you have only installed GLPK, type:

    ?> ./build-all.sh -DGLPK_PATH=<path in the system>

    Similarly, if you only want to use XPRESS, you don't need to put
information about GLPK.

    If you want to clean all folders, just run the shell script called
'clean-all.sh'.


=== AGSol installation on Ubuntu 14.04 ====
    
    As an example of how to install the AGSol solution in your machine, we
present the installation procedure in a clean Ubuntu 14.04.

(1) Install all the auxiliary dependencies, using the following command 
line:
?> sudo apt-get install cmake g++ libgmp-dev libmpfr-dev libboost-dev libboost-thread-dev

(2) Download the CGAL package using the following link:
http://gforge.inria.fr/frs/download.php/file/29125/CGAL-3.9.tar.gz.

(3) Enter in 'CGAL-3.9' folder and use the following commands:
?> cmake CMakeLists.txt
?> make
?> sudo make install

(4) Download the GLPK package from the link presented below:
http://ftp.gnu.org/gnu/glpk/glpk-4.52.1.tar.gz

(5) Install the AGSol package using the following commands:
?> ./configure
?> make
?> sudo make install

(6) Run the following command to configure the dynamic linker run-time 
bindings:
?> sudo ldconfig 

(7) Enter the agsol folder and compile AGSol running the command below: 
?> ./build-all.sh -DGLPK_PATH=/usr/local

The final binary will be placed in "agsol/art-gallery-pg" and is called 
'artGallerySolver'.

-------------------------------------------------------------------------------


4) How to run

    After compiling, enter the folder 'art-gallery-pg'. If the compilation was
successful, you will be able to see a binary file called 'artGallerySolver'.

    To execute, type the command below:

    ?> ./artGallerySolver <file .pol> <log file> <witness discretization> <solver mode>
    where,
        <file .pol> is the file representing the polygon. See the file format 
         in Section 4.1.  
        <log file> is the file where the summary log will be outputted. See the
         file format in Section 4.2.
        <witness discretization> consists in the technique chosen to select the
         initial Witness Set of our method. Current possibilities are:
            "ALL_VERTICES": includes all vertices of the polygon of input
            "CONVEX_VERTICES": includes all convex vertices of the polygon of
             input 
            "CHWA_POINTS": uses the method presented in [12] as a heuristic. The
             method assembles the initial witness set for our algorithm from
             the midpoints of all reflex-reflex edges and all convex vertices 
             from convex-reflex edges
            "CHWA_POINTS_EXTENDED": uses a combination of points in CHWA_POINTS 
             with the reflex vertices of the polygon
        <solver mode> is a positive integer that represents the mode used to
         solve SCP instances. Current possibilities are:
            "1": GLPK + Lagrangian Heuristic
            "2": GLPK
            "3": XPRESS + Lagrangian Heuristic
            "4": XPRESS

    If a final optimal solution is found, a file containing all of its
coordinates will be created in the same path of the instance file, but with a
.sol extension. The file format can be found in Section 4.3.

[12] Kyung-Yong Chwa, Byung-Cheol Jo, Christian Knauer, Esther Moet, René van
Oostrum, and Chan-Su Shin. 2006. Guarding art galleries by guarding witnesses.
Intern. Journal of Computational Geometry And Applications 16, 02n03 (2006),
205–226.


4.a) File Format: Instances
    Each file consists of one line divided in two parts. The first part
represents the outer boundary of the polygon (or the polygon itself, if there
are no holes). This part includes an initial integer value that represents the
number of vertices of the boundary and a counterclockwise sequence of the
vertices. Each vertex is represented by its x and y coordinates each of which
is written as the quotient of two integers (int/int).  

    The second part represents the holes of the polygon. First we have the
number of holes. After this, each hole is represented in the exact same way
described above for the boundary: an integer representing the size of the hole
and then a list of coordinates.

   As an example, here is the representation of a square with a hole (triangle)
inside: 
 4   1/1 1/1   100/2 1/1   500/10 50/1   1/1 100/2  1   3   10/1 100/10
80/2 10/1   100/4 40/1
         ______________
        |              |
        |      /\      |
        |     /  \     |
        |    /    \    |
        |   /      \   |
        |  /________\  |
        |______________|


4.b) File Format: Results
    Each file consists of one csv table with the following values: 
        -idFile: name of the instance
        -polSize: number of vertices of the instance
        -guards: size of the final solution
        -iterations: number of iterations necessary to find the optimal solution
        -initWitnSize: size of the initial Witness Set. See Section 4.4 in our 
         paper.
        -finalWitnSize: size of the final Witness Set.
        -initCandSize: size of the initial Guard Candidate Set. See Section 4.2
         of our paper.
        -finalCandSize: size of the final Guard Candidate Set.
        -maxHorizontal: maximum number of iterations necessary to solve an AGPFC 
         instance using Algorithm 3 (Section 4.2) in our paper.
        -IPSolved: number of SCP instances treated during execution.
        -initDisTime: Time spent creating initial discretization.
        -insertDisTime: Time spent inserting witnesses in the arrangement.
        -selectCandTime: Time spent selecting guard candidates.
        -initSolverTime: Time spent creating ILP matrix and initializing solver.
        -solverTime: Time spent solving all AGPFC instances.
        -scpResolTime: Time spent with ILP Solvers (including Lagrangian 
         Heuristic).
        -newDisTime: Time spent searching for new witnesses.
        -totalTime: Total CPU time spent by the program.


4.c) File Format: Solution
    Each file consists in 2 lines. The first one has the name of the instance
file that corresponds to the solution. The second line contains an integer
representing the size of the solution followed by a sequence of coordinates
representing the optimal guards placement.


4.d) The solution's viewer software:
    In our web page, it is possible to download a software which receives as
input both the polygon and its optimal solution and presents them in a friendly
user interface.

See our web page in the following URL:
 http://www.ic.unicamp.br/~cid/Problem-instances/Art-Gallery/AGPPG/

-------------------------------------------------------------------------------


5) Examples: 
    
?> ./artGallerySolver ../instances/agp2009a-simplerand/randsimple-100-1.pol log.txt CHWA_POINTS 1
Result: Runs software on simple instance with 100 vertices, using CHWA_POINTS
    initial discretization and GLPK Solver with Lagrangian Heuristic.

?> ./artGallerySolver ../instances/gD-ortho-ortho/gD_ortho-ortho_120:500v-50h_4.pol log-2.txt CONVEX_VERTICES 4
Result: Runs software on orthogonal polygon with holes with 500 vertices,
    using all convex vertices as initial discretization and XPRESS Solver.

Obs.: All of our instances that were experimented and whose results are
presented in the article can be found in the folder called 'instances'.

