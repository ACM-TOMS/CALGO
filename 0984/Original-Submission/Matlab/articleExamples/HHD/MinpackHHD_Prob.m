function Prob=MinpackHHD_Prob(vers)
% MinpackHHD_Prob: Minpack elastic-plastic torsion (EPT) problem data.
%
%     Computes start point, upper/lower bounds and problem dependent sigma 
%     parameters for the Human Heart Dipole (HHD) problem from the 
%     MINPACK-2 collection.
%
% USE:
%       Prob=MinpackHHD_Prob(vers)
% where
%   vers   : governs which version, 1 to 5, of the problem is used [1].
%
% Returns
%   Prob.x_0        : starting point
%   Prob.x_L        : lower bounds
%   Prob.x_U        : upper bounds
%   Prob.user.{summx,summy,suma,sumb,sumc,cumd,sume,sumf} : sigma
%   parameters for the problem.

% AUTHOR: S.A.Forth & K. Lenton
% DATE: 29/07/05
% Copyright 2005-2009: S.A. Forth, Cranfield University 
% REVISIONS:
% DATE  WHO   WHAT
% 30/6/09 SAF uses n instead of x and modified header
% 19/7/11 SAF no need to enter n

% Original Fortran Header comments follow
%      character*(*) task, prob
%      integer n, ldfjac
%   originally
%      double precision x(n), fvec(n), fjac(ldfjac,n)
%   changed to
%      double precision x(n), fvec(n), fjac(n,n)
% **********
%
% Subroutine dhhdfj
%
% This subroutine computes the function and the Jacobian matrix of
% the human heart dipole problem.
%
% The subroutine statement is
%
%   subroutine dhhdfj(n,x,fvec,fjac,ldfjac,task,prob)
%
% where
%
%   n is an integer variable.
%     On entry n is the number of variables. n = 8.
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
%         'XU'    Set x to the upper bound xu.
%
%     On exit task is unchanged.
%
%   prob is a character*5 variable.
%     On entry prob specifies the version of the problem. The
%        experiment label is the same as in Dennis, Gay, and Vu.
%
%               prob             experiment
%               ----             ---------
%              'DHHD1'             791129
%              'DHHD2'             791226
%              'DHHD3'              0121a
%              'DHHD4'              0121b
%              'DHHD5'              0121c
%
%     On exit prob is unchanged.
%
% MINPACK-2 Project. November 1993.
% Argonne National Laboratory and University of Minnesota.
% Brett M. Averick.
%
% **********

n = 8;
if nargin<1
    vers=1;
end

% set problem version dependent constants and x0
switch vers
    case 1
        Prob.user.summx = 0.485d0;
        Prob.user.summy = -0.0019d0;
        Prob.user.suma = -0.0581d0;
        Prob.user.sumb = 0.015d0;
        Prob.user.sumc = 0.105d0;
        Prob.user.sumd = 0.0406d0;
        Prob.user.sume = 0.167d0;
        Prob.user.sumf = -0.399d0;
        
        Prob.x_0 = [0.299d0; 0.186d0; -0.0273d0; 0.0254d0; -0.474d0; 0.474d0;...
            -0.0892d0; 0.0892d0];
    case 2
        Prob.user.summx = -0.69d0;
        Prob.user.summy = -0.044d0;
        Prob.user.suma = -1.57d0;
        Prob.user.sumb = -1.31d0;
        Prob.user.sumc = -2.65d0;
        Prob.user.sumd = 2.0d0;
        Prob.user.sume = -12.6d0;
        Prob.user.sumf = 9.48d0;
        Prob.x_0 = [-0.3d0; -0.39d0; 0.3d0; -0.344d0; -1.2d0; 2.69d0; 1.59d0; -1.5d0];
  
    case 3
        Prob.user.summx = -0.816d0;
        Prob.user.summy = -0.017d0;
        Prob.user.suma = -1.826d0;
        Prob.user.sumb = -0.754d0;
        Prob.user.sumc = -4.839d0;
        Prob.user.sumd = -3.259d0;
        Prob.user.sume = -14.023d0;
        Prob.user.sumf = 15.467d0;
        Prob.x_0 = [-0.041d0; -0.775d0; 0.03d0; -.047d0; -2.565d0; 2.565d0;...
             -0.754d0; 0.754d0];
    case 4
        Prob.user.summx = -0.809d0;
        Prob.user.summy = -0.021d0;
        Prob.user.suma = -2.04d0;
        Prob.user.sumb = -0.614d0;
        Prob.user.sumc = -6.903d0;
        Prob.user.sumd = -2.934d0;
        Prob.user.sume = -26.328d0;
        Prob.user.sumf = 18.639d0;
        Prob.x_0 = [-0.056d0; -0.753d0; 0.026d0; -0.047d0; -2.991d0; 2.991d0;...
             -0.568d0; 0.568d0];
    case 5
        Prob.user.summx = -0.807d0;
        Prob.user.summy = -0.021d0;
        Prob.user.suma = -2.379d0;
        Prob.user.sumb = -0.364d0;
        Prob.user.sumc = -10.541d0;
        Prob.user.sumd = -1.961d0;
        Prob.user.sume = -51.551d0;
        Prob.user.sumf = 21.053d0; 
        Prob.x_0 = [-0.074d0; -0.733d0; 0.013d0; -0.034d0; -3.632d0; 3.632d0;...
             -0.289d0; 0.289d0];
    case default
        error(['Illegal problem version vers=',num2str(vers)])
end
% x_L, x_U
Prob.x_L = -20*ones(n,1);
Prob.x_U = 20*[0; ones(n-1,1)];