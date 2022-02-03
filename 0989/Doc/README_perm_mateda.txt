
perm_mateda is a Matlab library for the solution of optimization problems defined on the permutation representation using estimation of distribution algorithms (EDAs). It implements probabilistic models specifically defined in the space of permutations. The library is conceived as an independent module that uses the Mateda framework, a general framework for solving optimization problems using EDAs.  

Getting started
----------------

0.- The first time you use perm_mateda, you need to download it (perm_mateda.zip) and unzip the file in the unit/folder of your preference.

Each time you want to use perm_mateda:

1.- Run MATLAB
2.- Inside MATLAB, change to the directory containing perm_mateda. It can be done by using the left panel (Current Folder), or throughout the command window (cd command). 
3.- run InitEnvironments.m by typing InitEnvironments in MATLAB command window)

Now you are ready to use perm_mateda. Different examples are provided, described in the section A quick introduction to perm_mateda. Three examples.

Main Features of perm_mateda
-----------------------------

- Implementation of Mallows and Generalized Mallows models [2-3].
- Implementation of Cayley, Kendall, and Ulam distances between permutations[1-5].
- Implementation of a number of optimization problems: Traveling Salesman Problem (TSP), Permutation Flowshop Scheduling Problem (PFSP), Linear Ordering Problem (LOP), and Quadratic Assignment Problem (QAP).
- As a control, the implementation of edge-histogram-based (EHM) and node-histogram-based (NHM) approaches to permutation problems [6] have been included.


Organization of the permutation directory
-----------------------------------------
- Consensus: Comprises the implementation of different methods for computing the consensus permutation given a set of permutations.            
- Distances: Implementation of Cayley, Kendall, and Ulam distances.         
- Histogram_Models: Contains the implementation of the EHM and NHM models. These are histogram-based models included for the sake of comparison with previous approaches.
- Mallows: Contains the implementation of the learning and sampling methods based on Mallows probabilistic models that uses different distances.       
- Operations: Contains a number of auxiliary functions, including two dependencies used for the generation of Ferrer Shapes. Programs colex.m and partition.m
- Problems: Implementation of TSP, PFSP, LOP, and QAP problems. It also contains test instances for these problems. 
- Scripts_Perm_Mateda: Contains three examples of Mallows EDAs using different parameters and applied on three different problems. It also contains post-processing steps for extracting and visualizing the results of the algorithms.   



A quick introduction to perm_mateda. Three examples
----------------------------------------------------

In the following we list the EDA running examples provided with perm_mateda. These examples are the scripts that have been used in the experimentation published in the paper. The scripts are included in the directory permutations/Scripts_Perm_Mateda. In order to run them, please move previously to the Scripts_Perm_Mateda directory, and type the name of the script as specified below (without extension) in MATLAB's command windows.

PLEASE NOTE that the examples have been configured to perform 30 repetitions and 500 generations per repetition. As this will take some time, for testing purposes the user can edit the example files, setting, for example, nrepetitions=1 and ngen=10.

- Example 1:  Application of Mallows EDA under the Ulam distance to the Linear Ordering Problem. The name of the file is Example_1_Mallows_Ulam_LOP.m
- Example 2:  Application of Generalized Mallows EDA under the Kendall distance to the Permutation Flowshop Scheduling Problem. The name of the file is Example_2_GMallows_Kendall_PFSP.m
- Example 3:  Application of Generalized Mallows EDA under the Cayley distance to the Quadratic Assignment Problem. The name of the file is Example_3_GMallows_Cayley_QAP.m



Some useful references
----------------------------------------------------------


[1] E. Irurozki, J. Ceberio, B. Calvo, J.A. Lozano. Mallows model under the Ulam distance: a feasible combinatorial approach. Neural Information Processing Systems 2014, Workshop on Analysis of Rank Data, Montreal, Canada 8-13 December 2014.
[2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011.
[3] J. Ceberio, E. Irurozki, A. Mendiburu, J.A. Lozano. A Distance-based Ranking Model Estimation of Distribution Algorithm for the Flowshop Scheduling Problem. IEEE Transactions on Evolutionary Computation. Vol 18, No. 2, Pp. 286-300. April 2014.
[4] J. Ceberio, E. Irurozki, A. Mendiburu, J.A. Lozano. A Review of Distances for the Mallows and Generalized Mallows Estimation of Distribution Algorithms. Journal of Computational Optimization and Applications. Vol. 62, No. 2, Pp. 545-564.
[5] J. Ceberio, R. Santana, A. Mendiburu, J.A. Lozano. Mixtures of Generalized Mallows models for solving the Quadratic Assignment Problem. 2015 IEEE Congress on Evolutionary Computation (CEC-2015),pp.2050-2057, Sendai, Japan, 25-28 May 2015.
[6] S. Tsutsui. Node histogram vs. edge histogram: A comparison of probabilistic model-building genetic algorithms in permutation domains. In: Evolutionary Computation, 2006. CEC 2006. IEEE Congress on. IEEE, 2006. p. 1939-1946.


