% NAVIER_FLOW
%
% Files
%   fpsetup_q0          - Q0 pressure convection-diffusion matrix 
%   fpsetup_q1          - Q1 pressure convection-diffusion matrix 
%   fpsetup_q2p1        - P1 pressure convection-diffusion matrix 
%   helpme_navier       - Navier-Stokes flow problem interactive help 
%   localbc_xycd        - imposes vector BC for Poisson error estimator 
%   navier_q1           - Q1 convection matrix 
%   navier_q2           - Q2 convection matrix 
%   navierpost_q1p0_bc  - postprocesses Poisson error estimator 
%   navierpost_q1p0_p   - computes Poisson error estimator for Q1-P0 
%   newton_q1           - Q1 convection derivative matrices 
%   newton_q2           - Q2 convection derivative matrices 
%   newtonbc            - imposes Dirichlet bc on Jacobian
%   null_pressure_index - index associated with the pressure nullspace
%   pressurebc          - fixes singularity in Laplacian matrix
%   solve_navier        - solve Navier-Stokes problem in square domain
%   solve_plate_navier  - solve Navier-Stokes problem in slit domain
%   solve_step_navier   - solve Navier-Stokes problem in step domain
%   Cpre_q1p0           - generate stabilizations for least sqrs commutator for Q1-P0 
%   Cpre_q1q1           - generate stabilizations for least sqrs commutator for Q1-Q1 
