% R_MOD Modified recurrence coefficients.
%
%    Given a weight function w through the first Nmax recurrence 
%    coefficients ab of its orthogonal polynomials, [abmod,Ncap,
%    kount]=R_MOD(N,ab) generates the first N recurrence
%    coefficients abmod of the modified weight function w/p, where
%    p is the polynomial of degree M defined by Eq.(3.1.75) 
%    incorporating the M poles selected. Their negative reciprocals
%    (the quantities zeta_mu in Eq.(3.1.71)) must be provided
%    in the first column of the global Mx2 array Z, the second 
%    column containing the respective multiplicities. It is assumed 
%    that complex poles occur in conjugate complex pairs and that
%    both poles of each pair are included in the array Z. Other 
%    global parameters are mc, mp, iq,idelta, irout, AB, eps0, Nmax, 
%    and M, which are used in the routine MCDIS called by R_MOD, and 
%    have the meaning explained there. Be sure that the variables 
%    ab, M, Z are declared global variables in the quadrature routine 
%    QUADRAT. The output variables Ncap and kount are inherited from 
%    the routine MCDIS with meanings explained there. The alpha- and
%    beta-coefficients of the given weight function are contained in
%    the first and second column of the Nmax x 2 input array ab, 
%    those of the modified weight function in the first and second
%    column of the Nx2 output array abmod.
%
%    See also MCDIS.
%
function [abmod,Ncap,kount]=r_mod(N,ab)
global mc mp iq idelta irout AB Z eps0 Nmax M
global theta c
if size(ab,1)<Nmax, error('input array ab too short'), end
if M>0
  [abmod,Ncap,kount]=mcdis(N,eps0,@quadrat,Nmax);
else
  Ncap=0; kount=0; abmod(1:N,:)=ab(1:N,:);
end
