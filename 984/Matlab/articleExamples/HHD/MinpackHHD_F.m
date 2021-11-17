function F=MinpackHHD_F(x,Prob)
% MinpackHHD_F: Computes function for the Human Heart Dipole (HHD) problem
%               from the MINPACK-2 collection.  
% USE:
%         F=MinpackHHD_F(x,Prob)
% where
%   x(8,1) : is an arbitrary vector x
%   Prob   : structure returned by MinpackHHD_Prob
%
% Returns
%   F        : function value

% AUTHOR: S.A.Forth & K. Lenton
% DATE: 26/06/06
% Copyright 2006-2009: S.A. Forth, Cranfield University
% REVISIONS:
% DATE  WHO   WHAT
% 30/6/09 SAF tidied some coding

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
n=length(x);
if n~=8
    error(['length(x) =',n,' is illegal, must = 8'])
end
% unpack constants
summx=Prob.user.summx ;
summy=Prob.user.summy ;
suma=Prob.user.suma ;
sumb=Prob.user.sumb ;
sumc=Prob.user.sumc ;
sumd=Prob.user.sumd ;
sume=Prob.user.sume ;
sumf=Prob.user.sumf ;
% function and Jacobian value
a = x(1);
b = x(2);
c = x(3);
d = x(4);
t = x(5);
u = x(6);
v = x(7);
w = x(8);
tv = t*v;
tt = t*t;
vv = v*v;
tsvs = tt - vv;
ts3vs = tt - 3*vv;
vs3ts = vv - 3*tt;
uw = u*w;
uu = u*u;
ww = w*w;
usws = uu - ww;
us3ws = uu - 3*ww;
ws3us = ww - 3*uu;
% % Evaluate the function
F=[a + b - summx;
    c + d - summy;
    t*a + u*b - v*c - w*d - suma;
    v*a + w*b + t*c + u*d - sumb;
    a*tsvs - 2*c*tv + b*usws - 2*d*uw - sumc;
    c*tsvs + 2*a*tv + d*usws + 2*b*uw - sumd;
    a*t*ts3vs + c*v*vs3ts + b*u*us3ws + d*w*ws3us - sume;
    c*t*ts3vs - a*v*vs3ts + d*u*us3ws - b*w*ws3us - sumf ];