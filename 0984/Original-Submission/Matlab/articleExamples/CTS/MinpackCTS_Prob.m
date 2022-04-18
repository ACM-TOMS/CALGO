function [Prob,nuse]=MinpackCTS_Prob()
% MinpackCTS_Prob: Minpack Coating Thickness Standardisation (CTS) problem data.
%
%    Sets data and scaling values for the Coating Thickness Standardisation 
%    (CTS) problem from the MINPACK-2 collection.
%
% USE:
%       [Prob,nuse] = MinpackCTS_Prob()
%
% Returns
%   Prob.x_0     : starting point
%   Prob.user.n  : number of variables to be fit
%   Prob.user.indvar : independent variables in fitting problem
%   Prob.user.y : dependent variables in fitting problem
%   Prob.user.scale1 : scaling for 
%   Prob.user.scale2 : scaling for 
%   nuse         : nuse=n
%
% See also MinpackCTS_F, MinpackCTS_FJ,MinpackCTS_Fnovec, MinpackCTS_FJnovec

% AUTHOR: S.A.Forth & K. Lenton
% DATE: 21/06/11
% Copyright 2011-2011: S.A. Forth, Cranfield University 
% REVISIONS:
% DATE  WHO   WHAT

% Original Fortran Header Comments follow  
%
%      subroutine dctsfj(m,n,x,fvec,fjac,ldfjac,task)
%      character*(*) task
%      integer m, n, ldfjac
%      double precision x(n), fvec(m), fjac(ldfjac,n)
% **********
%
% Subroutine dctsfj
%
% This subroutine computes the function and the Jacobian matrix of
% the coating thickness standardization problem.
%
% The subroutine statement is
%
%   subroutine dctsfj(m,n,x,fvec,fjac,ldfjac,task)
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
%
%     On exit task is unchanged.
%
% MINPACK-2 Project. November 1993.
% Argonne National Laboratory and University of Minnesota.
% Brett M. Averick.
%
% **********

% Initialization.
n = 134;

% Standard starting value
x_0 = zeros(n,1);
x_0(1:8) = [-8.0d0; 13.0d0; 1.2d0; 0.2d0; 0.1d0; 6.0d0; 5.5d0; -5.2d0];

% independent variables for data-fitting
indvar=[...
    0.714, 5.145; 0.7169, 5.241; 0.7232, 5.389; 0.7151, 5.211;
        0.6848, 5.154; 0.707, 5.105; 0.7177, 5.191; 0.7073, 5.013;
        0.6734, 5.582; 0.7174, 5.208; 0.7125, 5.142; 0.6947, 5.284;
        0.7121, 5.262; 0.7166, 6.838; 0.6894, 6.215; 0.6897, 6.817;
        0.7024, 6.889; 0.7026, 6.732; 0.68, 6.717; 0.6957, 6.468;
        0.6987, 6.776; 0.7111, 6.574; 0.7097, 6.465; 0.6809, 6.09;
        0.7139, 6.35; 0.7046, 4.255; 0.695, 4.154; 0.7032, 4.211;
        0.7019, 4.287; 0.6975, 4.104; 0.6955, 4.007; 0.7056, 4.261;
        0.6965, 4.15; 0.6848, 4.04; 0.6995, 4.155; 0.6105, 5.086;
        0.6027, 5.021; 0.6084, 5.04; 0.6081, 5.247; 0.6057, 5.125;
        0.6116, 5.136; 0.6052, 4.949; 0.6136, 5.253; 0.6032, 5.154;
        0.6081, 5.227; 0.6092, 5.12; 0.6122, 5.291; 0.6157, 5.294;
        0.6191, 5.304; 0.6169, 5.209; 0.5483, 5.384; 0.5371, 5.49;
        0.5576, 5.563; 0.5521, 5.532; 0.5495, 5.372; 0.5499, 5.423;
        0.4937, 7.237; 0.5092, 6.944; 0.5433, 6.957; 0.5018, 7.138;
        0.5363, 7.009; 0.4977, 7.074; 0.5296, 7.046];
    % y-values for data fitting
y=...
    [9.3636d0; 9.3512d0; 9.4891d0; 9.1888d0; 9.3161d0; 9.2585d0; 9.2913d0;... 
        9.3914d0; 9.4524d0; 9.4995d0; 9.4179d0; 9.468d0; 9.4799d0; 11.2917d0;...
        11.5062d0; 11.4579d0; 11.3977d0; 11.3688d0; 11.3897d0; 11.3104d0;...
        11.3882d0; 11.3629d0; 11.3149d0; 11.2474d0; 11.2507d0; 8.1678d0;...
        8.1017d0; 8.3506d0; 8.3651d0; 8.2994d0; 8.1514d0; 8.2229d0; 8.1027d0;...
        8.3785d0; 8.4118d0; 8.0955d0; 8.0613d0; 8.0979d0; 8.1364d0; 8.1700d0;...
        8.1684d0; 8.0885d0; 8.1839d0; 8.1478d0; 8.1827d0; 8.029d0; 8.1000d0;...
        8.2579d0; 8.2248d0; 8.2540d0; 6.8518d0; 6.8547d0; 6.8831d0; 6.9137d0;...
        6.8984d0; 6.8888d0; 8.5189d0; 8.5308d0; 8.5184d0; 8.5222d0; 8.5705d0;...
        8.5353d0; 8.5213d0; 8.3158d0; 8.1995d0; 8.2283d0; 8.1857d0; 8.2738d0;...
        8.2131d0; 8.2613d0; 8.2315d0; 8.2078d0; 8.2996d0; 8.3026d0; 8.0995d0;...
        8.2990d0; 9.6753d0; 9.6687d0; 9.5704d0; 9.5435d0; 9.6780d0; 9.7668d0;...
        9.7827d0; 9.7844d0; 9.7011d0; 9.8006d0; 9.7610d0; 9.7813d0; 7.3073d0;...
        7.2572d0; 7.4686d0; 7.3659d0; 7.3587d0; 7.3132d0; 7.3542d0; 7.2339d0;...
        7.4375d0; 7.4022d0; 10.7914d0; 10.6554d0; 10.7359d0; 10.7583d0;...
        10.7735d0; 10.7907d0; 10.6465d0; 10.6994d0; 10.7756d0; 10.7402d0;...
        10.6800d0; 10.7000d0; 10.8160d0; 10.6921d0; 10.8677d0; 12.3495d0;...
        12.4424d0; 12.4303d0; 12.5086d0; 12.4513d0; 12.4625d0; 16.2290d0;...
        16.2781d0; 16.2082d0; 16.2715d0; 16.2464d0; 16.1626d0; 16.1568d0];
    % scaling parameters
scale1=4.08d0;
scale2=0.417d0;


% store required values in Prob
Prob.x_0 = x_0;
Prob.user.n = n;
Prob.user.indvar = indvar;
Prob.user.y = y;
Prob.user.scale1 = scale1;
Prob.user.scale2 = scale2;
nuse=n;