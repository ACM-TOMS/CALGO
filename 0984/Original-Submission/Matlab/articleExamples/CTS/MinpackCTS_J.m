function J=MinpackCTS_J(x,Prob)
% MinpackCTS_J: Minpack Coating Thickness Standardisation (CTS) Jacobian.
%
%    Computes Jacobian for the Coating Thickness
%    Standardisation (CTS) problem from the MINPACK-2 collection.
%
% USE:
%           J = MinpackCTS_J(x,Prob)
% where
%   x    : solution vector with length(x)=134
%   Prob : structure created by MinpackEPT_Prob with components
%     Prob.user.n      : number of variables to be fit
%     Prob.user.indvar : independent variables in fitting problem
%     Prob.user.y      : dependent variables in fitting problem
%     Prob.user.scale1 : scaling for
%     Prob.user.scale2 : scaling for

% AUTHOR: S.A.Forth & K. Lenton
% DATE: 21/06/11
% Copyright 2011-2011: S.A. Forth, Cranfield University
% REVISIONS:
% DATE  WHO   WHAT

% Original Fortran Header Comments follow
%
%      subroutine dctsfj(m,n,x,fvec,J,ldJ,task)
%      character*(*) task
%      integer m, n, ldJ
%      double precision x(n), fvec(m), J(ldJ,n)
% **********
%
% Subroutine dctsfj
%
% This subroutine computes the function and the Jacobian matrix of
% the coating thickness standardization problem.
%
% The subroutine statement is
%
%   subroutine dctsfj(m,n,x,fvec,J,ldJ,task)
%
% where
%
%   m is an integer variable.
%     On entry m is the number of functions.
%        For the coating thickness standardization problem m = 252.
%     On exit m is unchanged.
%
%   n is an integer variable.
%     On entry n is the number of variables.
%        For the coating thickness standardization problem n = 134.
%     On exit n is unchanged.
%
%   x is a double precision array of dimension n.
%     On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
%        Otherwise x need not be specified.
%     On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
%        x is set according to task.
%
%   fvec is a double precision array of dimension m.
%     On entry fvec need not be specified.
%     On exit fvec contains the function evaluated at x if
%        task = 'F' or 'FJ'.
%
%   J is a double precision array of dimension (ldJ,n).
%     On entry J need not be specified.
%     On exit J contains the Jacobian matrix evaluated at x if
%        task = 'J' or 'FJ'.
%
%   ldJ is an integer variable.
%      On entry ldJ is the leading dimension of J.
%      On exit ldJ is unchanged.
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
%
%     On exit task is unchanged.
%
% MINPACK-2 Project. November 1993.
% Argonne National Laboratory and University of Minnesota.
% Brett M. Averick.
%
% **********

if nargin~=2
    error ('MADMinpack:CTS:MinpackCTS_J:nargin',...
        'MADMinpack:CTS:MinpackCTS_J:nargin - must have 2 input arguments x and Prob')
end

n=length(x);
if (n ~= Prob.user.n)
    error ('MADMinpack:CTS:MinpackCTS_J:n',...
        ['MADMinpack:CTS:MinpackCTS_J:n - length(x) must be ',num2str(Prob.user.n)])
end

% number of equations
m = 2*(n-8);
mdiv4=m/4;

% unpack user parameters
indvar = Prob.user.indvar;
scale1 = Prob.user.scale1;
scale2 = Prob.user.scale2;

J = sparse(m,n);

J(1:mdiv4,1:4) = [ones(mdiv4,1), indvar(:,1) + x(9:8+mdiv4), indvar(:,2) + x(mdiv4+9:8+2*mdiv4),...
    (indvar(:,1)+x(9:8+mdiv4)).*(indvar(:,2)+x(mdiv4+9:8+2*mdiv4))];
J(1:mdiv4,9:8+mdiv4) = spdiags(x(2) + x(4)*(indvar(:,2)+x(mdiv4+9:8+2*mdiv4)),0,sparse(mdiv4,mdiv4));
J(1:mdiv4,9+mdiv4:8+2*mdiv4) = spdiags(x(3) + x(4)*(indvar(:,1)+x(9:8+mdiv4)),0,sparse(mdiv4,mdiv4));

J(mdiv4+1:2*mdiv4,5:8) = [ones(mdiv4,1), indvar(:,1) + x(9:8+mdiv4), indvar(:,2) + x(mdiv4+9:8+2*mdiv4),...
    (indvar(:,1)+x(9:8+mdiv4)).*(indvar(:,2)+x(mdiv4+9:8+2*mdiv4))];
J(mdiv4+1:2*mdiv4,9:8+mdiv4) = spdiags(x(6) + x(8)*(indvar(:,2)+x(mdiv4+9:8+2*mdiv4)),0,sparse(mdiv4,mdiv4));
J(mdiv4+1:2*mdiv4,9+mdiv4:8+2*mdiv4) = spdiags(x(7) + x(8)*(indvar(:,1)+x(9:8+mdiv4)),0,sparse(mdiv4,mdiv4));
J(2*mdiv4+1:3*mdiv4,9:8+mdiv4) = spdiags(scale1(ones(mdiv4,1)),0,sparse(mdiv4,mdiv4));
J(3*mdiv4+1:4*mdiv4,9+mdiv4:8+2*mdiv4) = spdiags(scale2(ones(mdiv4,1)),0,sparse(mdiv4,mdiv4));
