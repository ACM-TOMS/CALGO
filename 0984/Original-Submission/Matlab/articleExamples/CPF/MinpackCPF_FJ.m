function [F,J]=MinpackCPF_FJ(x,Prob)
% MinpackCPF_FJ: Function & Jacobian for Minpack CPF problem
%
%     Computes function and, optionally, the Jacobian for the combustion of
%     propane full formulation (CPF) problem from the MINPACK-2 collection.
%
% USE:
%       F=MinpackCPF_FJ(x)
% where
%   x      : input vector with length(x)=11
%
% Returns
%  F  : Function
%  J  : Jacobian

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
if nargin~=2
    error ('MADMinpack:CPF:MinpackCPF_FJ:nargin',...
        'MADMinpack:CPF:MinpackCPF_FJ:nargin - must have 2 inputs')
end

% check n
n=length(x);
if n ~= 11
    error ('MADMinpack:CPF:MinpackCPF_Prob:n',...
        ['MADMinpack:CPF:MinpackCPF_Prob:n - n = length(x) = ',num2str(n),' is illegal, n must = 11'])
end

% constants
k=Prob.user.k;
p=Prob.user.p;
rr=Prob.user.rr;
sqrtp = sqrt(p);

% coding
pdx = p/x(11);
sqpdx = sqrt(pdx);
xtau = sum(x(1:n-1));

F = [x(1) + x(4) - 3.0d0;...
    2*x(1) + x(2) + x(4) + x(7) + x(8) + x(9) + 2*x(10) - rr;...
    2*x(2) + 2*x(5) + x(6) + x(7) - 8.0d0;...
    2*x(3) + x(9) - 4.0d0*rr;...
    k(5)*x(2)*x(4) - x(1)*x(5);...
    k(6)*sqrt(x(2)*x(4)) - sqrt(x(1))*x(6)*sqpdx;...
    k(7)*sqrt(x(1)*x(2)) - sqrt(x(4))*x(7)*sqpdx;...
    k(8)*x(1) - x(4)*x(8)*pdx;...
    k(9)*x(1)*sqrt(x(3)) - x(4)*x(9)*sqpdx;...
    k(10)*x(1)^2 - (x(4)^2)*x(10)*pdx;
    x(11) - xtau];

if nargout==2
    J=zeros(n);
    J(n,1:n-1)=-1;
    J(n,n) = 1;
    xfrac = 1/(sqrt(x(11))^3);
    
    J(1,1) = 1;
    J(2,1) = 2;
    J(5:10,1) = [-x(5); -0.5*x(6)*sqpdx/sqrt(x(1)); 0.5*k(7)*sqrt(x(2))/sqrt(x(1)); k(8); k(9)*sqrt(x(3)); 2*k(10)*x(1)];
    
    J(2,2) = 1;
    J(3,2) = 2;
    J(5:7,2) = [k(5)*x(4); 0.5*k(6)*sqrt(x(4))/sqrt(x(2)); 0.5*k(7)*sqrt(x(1))/sqrt(x(2))];
    
    J(4,3) = 2;
    J(9,3) = 0.5*k(9)*x(1)/sqrt(x(3));
    
    J(1:2,4) = 1;
    J(5:10,4) = [k(5)*x(2); 0.5*k(6)*sqrt(x(2))/sqrt(x(4)); -0.5*x(7)*sqpdx/sqrt(x(4)); -x(8)*pdx; -x(9)*sqpdx; -2*x(4)*x(10)*pdx];
    
    J(3,5) = 2;
    J(5,5) = -x(1);
    
    J(3,6) = 1;
    J(6,6) = -sqrt(x(1))*sqpdx;
    
    J(2:3,7) = 1;
    J(7,7) = -sqrt(x(4))*sqpdx;
    
    J(2,8) = 1;
    J(8,8) = -x(4)*pdx;
    
    J(2,9) = 1;
    J(4,9) = 1;
    J(9,9) = -x(4)*sqpdx;
    
    J(2,10) = 2;
    J(10,10) = -(x(4)^2)*pdx;
    
    J(6:10,11) = [0.5*sqrt(x(1))*x(6)*sqrtp*xfrac; 0.5*sqrt(x(4))*x(7)*sqrtp*xfrac; x(4)*x(8)*p/(x(11)^2); 0.5*x(4)*x(9)*sqrtp*xfrac; x(4)^2*x(10)*p/(x(11)^2)];
    
end