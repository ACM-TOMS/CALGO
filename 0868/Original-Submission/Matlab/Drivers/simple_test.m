%simple_test
%This is a matlab script where we will illustrate how to
%use the three codes: da0glob.m, da2glob.m, da3glob.m and in addition quadl.
%Just write simple_test in a Matlab's command window to run this code.

%   Terje O. Espelid 31/01/06

% Give the interval

format long
 a=0;b=1;

% and the tolerance
 tol=10^-12;

% compute the integral to an ABSOLUTE error tolerance tol using da0glob
 [q_estimate(1),nfun(1)]=da0glob(@simple_example,[a,b],tol);

% compute the integral to a relative error tolerance to using da2glob.
 [q_estimate(2),nfun(2)]=da2glob(@simple_example,[a,b],0,tol);

% compute the integral to an absolute error tolerance and apply trace using
% da3glob.
 [q_estimate(3),nfun(3)]=da3glob(@simple_example,[a,b],tol,0,1);

% compute the integral to an absolute error tolerance using quadl
 [q_estimate(4),nfun(4)]=quadl(@simple_example,a,b,tol);

% The exact answer to this problem
 exact=.013492485649467772692
% The four results:
 q_estimate
% The number of function evaluations needed.
 nfun
%to achieve the tolerance
tol
% The error in the four cases
 error_q=q_estimate-exact

% This gives the following output (Dell Computer OPTIPLEX GX260)
% First the trace results from da3glob.m: 21 different intervals, the code
% has used all rules in the process
% but for the final intervals the rule has used the 9, 17 and 33 point
% rules only. Totally this code has used 337 evaluation points.
%
%    left_endpoint    right_endpoint     quadrature value   rule used
%-------------------------------------------------------------------------
%Here is the output: first the trace output from da3glob; then the results
%>>simple_test
%                  0   0.06250000000000   0.00013300798657  17.00000000000000

%   0.50000000000000   0.75000000000000   0.00002063797601  17.00000000000000

%   0.25000000000000   0.37500000000000   0.00008074659032  17.00000000000000

%   0.12500000000000   0.12695312500000   0.00096016631264  17.00000000000000

%   0.18750000000000   0.21875000000000   0.00011675002319  17.00000000000000

%   0.15625000000000   0.17187500000000   0.00027095468190  17.00000000000000

%   0.14062500000000   0.14843750000000   0.00072308593441  17.00000000000000

%   0.06250000000000   0.09375000000000   0.00023502342298  17.00000000000000

%   0.13281250000000   0.13671875000000   0.00202213906880  17.00000000000000

%   0.09375000000000   0.10937500000000   0.00037227411896  17.00000000000000

%   0.10937500000000   0.11718750000000   0.00049363836407  17.00000000000000

%   0.12890625000000   0.13085937500000   0.00189310302750  17.00000000000000

%   0.11718750000000   0.12500000000000   0.00155483848283  33.00000000000000

%   0.37500000000000   0.50000000000000   0.00002613800912  17.00000000000000

%   0.13671875000000   0.14062500000000   0.00087818280915  17.00000000000000

%   0.14843750000000   0.15625000000000   0.00030485459307  17.00000000000000

%   0.17187500000000   0.18750000000000   0.00012387836645  17.00000000000000

%   0.13085937500000   0.13281250000000   0.00175257930287  17.00000000000000

%   0.75000000000000   1.00000000000000   0.00000877161015   9.00000000000000

%   0.21875000000000   0.25000000000000   0.00005584101273   9.00000000000000

%   0.12695312500000   0.12890625000000   0.00146587395577  17.00000000000000


%exact =

%   0.01349248564947


%q =

%   0.01349248564949   0.01349248564947   0.01349248564949   0.01349248564947


%nfun =

%   401   401   337   888


%tol =

%     1.000000000000000e-12


%error =

%   1.0e-13 *

%   0.20939847078516   0.00045102810375   0.21630267021955   0.00086736173799
