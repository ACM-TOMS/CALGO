function [Prob,nuse]=MinpackCPF_Prob()
% MinpackCPF_Prob: Minpack combustion of propane full formulation (CPF) problem data.
%
%     Computes start point, lower bounds and problem dependent parameters 
%     for the combustion of propane full formulation (CPF) problem from the
%     MINPACK-2 collection.
%
% USE:
%       [Prob,nuse]=MinpackCPF_Prob()
%
% Returns
%   Prob.x_0        : starting point
%   Prob.x_L        : lower bounds
%   Prob.x_U        : upper bounds
%   Prob.user.{summx,summy,suma,sumb,sumc,cumd,sume,sumf} : sigma
%   parameters for the problem.
%   nuse            : value of n used

% AUTHOR: S.A.Forth & K. Lenton
% DATE: 30/06/09
% Copyright 2009-2009: S.A. Forth, Cranfield University
% REVISIONS:
% DATE  WHO   WHAT

% Original Fortran Header Comments follow
%
% subroutine dcpffj(n,x,fvec,fjac,ldfjac,task)
%      character*(*) task
%      integer n, ldfjac
%      double precision x(n), fvec(n), fjac(ldfjac,n)
% **********
%
% Subroutine dcpffj
%
% This subroutine computes the function and the Jacobian matrix of
% the combustion of propane (full formulation) problem.
%
% The subroutine statement is
%
%   subroutine dcpffj(n,x,fvec,fjac,ldfjac,task)
%
% where
%
%   n is an integer variable.
%     On entry n is the number of variables. n = 11.
%     On exit n is unchanged.
%
%   x is a double precision array of dimension n.
%     On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
%        Otherwise x need not be specified.
%     On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
%        x is set according to task.
%
%   fvec is a double precision array of dimension n.
%     On entry fvec need not be specified.
%     On exit fvec contains the function evaluated at x if
%        task = 'F' or 'FJ'.
%
%   fjac is a double precision array of dimension (ldfjac,n).
%     On entry fjac need not be specified.
%     On exit fjac contains the Jacobian matrix evaluated at x if
%        task = 'J' or 'FJ'.
%
%   ldfjac is an integer variable.
%      On entry ldfjac is the leading dimension of fjac.
%      On exit ldfjac is unchanged.
%
%   task is a character variable.
%     On entry task specifies the action of the subroutine:
%
%        task               action
%        ----               ------
%         'F'     Evaluate the function at x.
%         'J'     Evaluate the Jacobian matrix at x.
%         'FJ'    Evaluate the function and the Jacobian at x.
%         'XS'    Set x to the standard starting point xs.
%         'XL'    Set x to the lower bound xl.
%
%     On exit task is unchanged.
%
% MINPACK-2 Project. November 1993.
% Argonne National Laboratory and University of Minnesota.
% Brett M. Averick.
%
% **********

% check we have correct number of arguments
if nargin>=1
    error ('MADMinpack:CPF:MinpackCPF_Prob:nargin',...
        'MADMinpack:CPF:MinpackCPF_Prob:nargin - must have less than one input')
end
n=11;

% check n
if n ~= 11
    error ('MADMinpack:CPF:MinpackCPF_Prob:n',...
        ['MADMinpack:CPF:MinpackCPF_Prob:n - n = ',num2str(n),' is illegal, n must = 11'])
end
nuse=11;

% constants
k=[0.0d0 0.0d0 0.0d0 0.0d0 1.930d-1 2.597d-3 3.448d-3 1.799d-5 2.155d-4 3.846d-5];
p=4.0d1;
sqrtp = sqrt(p);
Prob.user.k=k;
Prob.user.p=p;
Prob.user.rr=1.0d1;

% bounds and start point
Prob.x_L= zeros(n,1);
Prob.x_U = [];
Prob.x_0 = [5.0d0; 2.5d0; 5.0d0; 1.0d-1; 5.0d-2*k(5); k(6)/sqrtp; 5.0d1*k(7)/sqrtp;...
            1.0d3*k(8)/p; 5.0d2*k(9)/sqrtp; 5.0d4*k(10)/p; 2.0d1];