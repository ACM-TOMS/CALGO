% SOLVERS
%
% Files
%   a_cdt              - matrix-vector product for convection-diffusion operator
%   a_nst              - matrix-vector product for linearized Navier-Stokes operator
%   bicgstab_ell       - bicgstab(ell) for right-preconditioned iterates
%   bicgstab_ell_r     - bicgstab(ell) iteration with right preconditioning 
%   cg_test            - CG convergence demo 
%   givapp             - apply a sequence of Givens rotations
%   gmres_r            - GMRES iteration with right preconditioning 
%   helpme_it          - iterative solvers interactive help 
%   helpme_mg          - geometric multigrid interactive help
%   ilu0               - incomplete factorization with no fill of sparse matrix%   it_solve           - driver for iterative solution of predefined problem
%   m_amgt             - AMG preconditioner for convection-diffusion operator
%   m_bfbt             - ideal least squares commutator preconditioner (unscaled)
%   m_diagt            - action of diagonal preconditioning operator
%   m_fp               - ideal pressure convection-diffusion preconditioner
%   m_fp_amg           - AMG iterated pressure convection-diffusion preconditioner
%   m_fp_mg            - GMG iterated pressure convection-diffusion preconditioner
%   m_ilut             - incomplete LU preconditioner
%   m_mg               - GMG preconditioner for scalar problems
%   m_nonet            - "no preconditioning" operator
%   m_st_block         - block preconditioner for Stokes equations
%   m_st_mg            - block MG preconditioner for Stokes equations
%   m_sxbfbt           - ideal stabilized least squares commutator preconditioner
%   m_xbfbt            - ideal least squares commutator preconditioner
%   m_xbfbt_amg        - AMG iterated least squares commutator preconditioner
%   m_xbfbt_mg         - GMG iterated least squares commutator preconditioner
%   mg_apsetup_q1      - Q1 pressure diffusion matrix for GMG
%   mg_cavity_domain   - square cavity Q2 grid generator for GMG
%   mg_cd              - preconditioner for convection-diffusion problem
%   mg_cd_setup        - GMG convection-diffusion problem on square domain
%   mg_diff            - GMG preconditioner for diffusion problem
%   mg_diff_q2q1grid   - Q2-Q1 element grid generator for GMG
%   mg_diff_setup      - GMG diffusion problem on square domain
%   mg_diff_setup_ell  - GMG diffusion problem on L-shaped domain
%   mg_diff_setup_step - GMG diffusion problem on step domain
%   mg_ellblock        - prolongation for part of L-shaped and step domains
%   mg_iter            - performs one GMG iteration
%   mg_ns              - GMG preconditioner for Navier-Stokes equations
%   mg_ns_cd           - GMG for convection-diffusion problem  (Navier-Stokes)
%   mg_ns_cd_setup     - GMG convection-diffusion problem (Navier-Stokes)
%   mg_ns_diff         - GMG for diffusion problem (Navier-Stokes)
%   mg_ns_diff_lstsq   - GMG for diffusion problem (Navier-Stokes, least-squares)
%   mg_ns_iter         - performs one MG iteration (Navier-Stokes)
%   mg_ns_prolong_step - GMG prolongation operator for step domain (Navier-Stokes)
%   mg_ns_q1cd         - convection-diffusion matrix generator for GMG (Navier-Stokes)
%   mg_ns_q1cd_supg    - streamline diffusion matrix generator for GMG (Navier-Stokes)
%   mg_ns_smooth       - smoothers for GMG iteration (Navier-Stokes)
%   mg_post            - postsmoothing for GMG
%   mg_pre             - presmoothing for GMG
%   mg_prolong         - GMG prolongation operator for square domain
%   mg_prolong_ell     - GMG prolongation operator for L-shaped domain
%   mg_prolong_step    - GMG prolongation operator for step domain
%   mg_q1cd            - convection-diffusion matrix generator for GMG
%   mg_q1cd_supg       - streamline diffusion matrix generator for GMG 
%   mg_q1diff          - bilinear diffusion matrix generator for GMG 
%   mg_q1grid          - bilinear element grid generator for GMG
%   mg_smooth          - smoothers for GMG on square domain
%   mg_smooth_ell      - smoothers for GMG on L-shaped domain
%   mg_smooth_step     - smoothers for GMG on step domain
%   mg_solve           - driver for GMG solution of predefined problem 
%   mg_step_domain     - Q2 grid generator for GMG on step domain
%   mg_zerobc          - imposes zero boundary conditions
%   resplot            - plot residuals computed by iterative solvers
