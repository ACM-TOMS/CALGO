% STOKES_FLOW
%
% Files
%   flowbc             - imposes inflow boundary condition 
%   forwardstep_stokes - set up flow problem in forward step domain 
%   helpme_stokes      - Stokes flow problem interactive help 
%   infsup             - computes inf-sup eigenvalue distribution 
%   localbc_xy         - imposes vector BC for Poisson error estimator 
%   longstep_stokes    - set up flow problem in extended step domain 
%   pipe_stokes        - set up inflow/outflow problem in square domain 
%   plate_stokes       - set up flow problem in slit domain 
%   q1div              - computes norm of divergence of Q1 flow solution
%   q2div              - computes norm of divergence of Q2 flow solution
%   quad_stokes        - set up Stokes problem in quadrilateral domain 
%   solve_step_stokes  - solve Stokes problem in step domain
%   solve_stokes       - solve Stokes problem in square domain
%   specific_flow      - Reference problem 5.3 default inflow condition 
%   square_stokes      - set up flow problem in unit square domain 
%   step_stokes        - set up flow problem in standard step domain 
%   stokes_q1p0        - vectorized Q1-P0 matrix generator
%   stokes_q1q1        - vectorized Q1-Q1 matrix generator
%   stokes_q2p1        - vectorized Q2-P1 matrix generator
%   stokes_q2q1        - vectorized Q2-Q1 matrix generator
%   stokespost_q1p0_bc - postprocesses Poisson error estimator 
%   stokespost_q1p0_p  - computes Poisson error estimator for Q1-P0 
%   stokespost_q1q1_bc - postprocesses Poisson error estimator 
%   stokespost_q1q1_p  - computes Poisson error estimator for Q1-Q1 
%   stressjmps_q1p0    - stress jumps for rectangular Q1-P0 grid 
%   stressjmps_q1q1    - stress jumps for rectangular Q1-Q1 grid 
