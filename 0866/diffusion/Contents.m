% DIFFUSION
%
% Files
%   deriv        - evaluates derivatives of bilinear shape functions 
%   diffpost_bc  - postprocesses local Poisson error estimator 
%   diffpost_p   - computes local Poisson error estimator for Q1 solution 
%   diffpost_res - computes Q1 element residual error estimator 
%   ell_diff     - solve Poisson problem in L-shaped domain 
%   femq1_diff   - vectorized bilinear coefficient matrix generator
%   femq2_diff   - vectorized biquadratic coefficient matrix generator
%   gauss_source - evaluates source term at Gauss point 
%   helpme_diff  - diffusion problem interactive help 
%   lderiv       - evaluates derivatives of linear shape functions 
%   localbc_p    - imposes Dirichlet BC for Poisson error estimator 
%   lshape       - evaluates linear shape functions 
%   nonzerobc    - imposes Dirichlet boundary condition 
%   q1fluxjmps   - computes flux jumps for rectangular Q1 grid 
%   q1res_diff   - computes interior residuals for rectangular Q1 grid
%   qderiv       - evaluates derivatives of biquadratic shape functions 
%   qshape       - evaluates biquadratic shape functions 
%   quad_diff    - solve Poisson problem in quadrilateral domain 
%   shape        - evaluates bilinear shape functions 
%   specific_bc  - (current) problem boundary condition 
%   specific_rhs - (current) problem forcing function
%   square_diff  - solve Poisson problem in unit square domain 