% Contents.m in JSR_louvain
%
% NOTE :  Functions with an * uses solve_semi_definite_program. 
%         By default, the solver is SeDuMi, but it is possible to use a 
%         different solver. See the help for more information.
%
% You can download SeDuMi here :
% http://perso.uclouvain.be/raphael.jungers/sites/default/files/sedumi.zip
%
% Upper-level methods
%
%    jsr                   - * Algorithm doing pre-processing and launching
%                            different methods to efficiently find
%                            bounds on the jsr
%    itMeth                - Computes all products of certain lengths and
%                            launches user-specified methods on those sets
%
% Demos
%
%    demo1_JSR             - An introduction on how to call functions,
%                            the kind of outputs and how to specify 
%                            values for the parameters
%    demo2_JSR             - * Explained model script for the launch of 
%                            two methods
%    demo3_JSR             - * Explains how to use jsr_pathcomplete method.
%
% Option generator
%
%    jsrsettings           - Creates a structure to be used with any
%                            method in order to specify parameters and 
%                            options
%
% Pre-processing routines
%
%    comp2real             - Represents a set of complex nxn matrices as
%                            a set of real (2*n)x(2*n) matrices
%    jointTriangul         - Tries heuristically to jointly block-triangularise 
%                            a set of matrices. Returns the sets of diagonal
%                            blocks and unitary change of basis matrix
%    permTriangul          - Looks for a permutation of the rows and
%                            columns of nonnegative matrices in a set
%                            that jointly block-triangularise them
%    quickElim             - Eliminates irrelevant blocks from a set of
%                            diagonal blocks
%
%                               
% Methods 
%
%   Conic
%    jsr_conic_ellipsoid                - * Approximates the jsr using
%                                         ellipsoidal norms
%    jsr_conic_linear                   - * Approximates the jsr of a set of  
%                                         nonnegative matrices using the 
%                                         joint conic radius in the 
%                                         positive orthant
%    jsr_pathcomplete                   - * Approximate the jsr using
%                                         piecewise quadratic Lyapunov
%                                         functions, with the help of a
%                                         pathcomplete graph.
%
%   Lift
%    jsr_lift_semidefinite              - Approximates the jsr using
%                                         semidefinite liftings  
%    
%   Norm
%    jsr_norm_balancedComplexPolytope   - Approximates the jsr using b.c.p.'s
%    jsr_norm_balancedRealPolytope      - Approximates the jsr using b.r.p.'s
%    jsr_norm_conitope                  - * Approximates the jsr using lifted
%                                         BCP
%    jsr_norm_linearRelaxation2D        - Approximates the jsr using
%                                         Linear-Relaxation in 2D
%                                         (Heuristic)
%    jsr_norm_maxRelaxation             - Approximates the jsr using
%                                         Max-Relaxation. (Heuristic)
%    jsr_norm_maxRelaxation2D           - Approximates the jsr using
%                                         Max-Relaxation in 2D (Heuristic)
%   
%   Optimization
%    jsr_opti_sos                       - * Approximates the jsr using sum of
%                                         squares
%    
%   Product
%    jsr_prod_bruteForce                - Approximates the jsr using brute force
%    jsr_prod_Gripenberg                - Approximates the jsr using branch
%                                         and bound
%    jsr_prod_lowerBruteForce           - Gives a lower bound on the jsr using brute force
%    jsr_prod_pruningAlgorithm          - Approximates the jsr using pruning algorithm
%
% Benchmark
%    BTVMat                - Generates set of matrices introduced by
%                            Blondel, Theys and Vladimirov. Related to finiteness 
%                            conjecture
%    overlapMat            - Generates set of matrices whose JSR
%                            characterizes the asymptotic growth rate of 
%                            overlap-free words  
%    waveletMat            - Generates set of matrices whose JSR is related  
%                            to the continuity of Daubechies' wavelets
%
% Subroutines
%    available_memory      - Estimates free memory on the computer
%    buildProduct          - Computes the product corresponding to a
%                            specified sequence of indices
%    cellDivide            - Divides each matrix in a cell array by a scalar
%    com_eig               - Computes a common eigenvector of two matrices
%    debruijn              - Builds a DeBruijn graph for a set of label
%    deperiod              - Finds a periodic sequence of indices
%    find_pathcomplete_lyapunov - * Solves a quazi-convex optimisation problem
%                            for jsr_pathcomplete.
%    findRow               - Finds a given row in an ordered matrix
%    genNecklaces          - Generation of all n-bead necklaces with k 
%                            colors
%    genPerms              - Generation of all non-negative integer
%                            n-tuples of given (max) sum
%    generate_graph        - Generates graphs needed for some methods
%    generate_pathcomplete_SDP - Generates a SDP problem associated to
%                            a given graph for find_pathcomplete_lyapunov
%    graphSCC              - Finds the strongly connected components of 
%                            a directed graph using Tarjan's algorithm
%    iscellcomplex         - Returns 1 if there is complex matrices in set M
%    liftProduct           - Generation of the set of all k-products 
%                            of matrices of a set
%    msg                   - Prints a message with arguments in command
%                            line and in a logFile
%    polyliftedNorm        - * Computes polytope norm of a symmetric SDP 
%                            matrix w.r.t. a given essential set
%    pruneSet              - Removes majorated matrices in a set of SDP
%                            matrices 
%    tens2cell             - Converts a tensor to a cell array of matrices
%    solve_semi_definite_program - * Solves a standard semi-definite-program
%                            with specific outputs
%
%