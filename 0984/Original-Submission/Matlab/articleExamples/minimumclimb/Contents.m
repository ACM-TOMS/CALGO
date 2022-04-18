% Minimum Time to Climb Optimal Control Problem
%
% Files:
% minimumClimbAuxdata.m      - builds auxiliary data structure for problem
% minimumClimbConstraints.m  - Constraint function for the NLP
% minimumClimbDynamics.m     - Vectorized dynamics file
% minimumClimbObjective.m    - Objective function for the NLP
% minimumClimbRatios.m       - Compares ADiGator/ADiMat/INTLAB/MAD ratios
% minimumClimbSolve_nonvec.m - Solves the optimal control problem using
%                              ADiGator and IPOPT, without vectorization
% minimumClimbSolve_vec.m    - Solves the optimal control problem using
%                              ADiGator and IPOPT, with vectorization
% minimumClimbSum.m          - Takes the sum of the dynamics file so that
%                              ADiMat can be using in FM over RM for second 
%                              derivative.
% setupMinimumClimb.m        - Sets up the LGR collocation, builds initial
%                              guess for NLP.
% a_minimumClimbSum.m        - ADiMat generated reverse mode file
% g_minimumClimbDynamics.m   - ADiMat gererated forward mode file
%                              (slightly modified to allow for cells)
%
% Directories:
% ipopt  - IPOPT help and IPOPT MEX files.