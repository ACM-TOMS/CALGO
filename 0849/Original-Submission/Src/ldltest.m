%
% ldltest.m:  compile and test the LDL mexFunction
%
% LDL Version 1.1 (Apr. 22, 2005), Copyright (c) 2005 by Timothy A Davis,
% University of Florida.  All Rights Reserved.  See README for the License.

clear
help ldl
mex -inline ldl.c ldlmex.c

A = sparse ([ ]) ;

[L, D, Parent, fl] = ldl (A) ;
[L, D, Parent, fl] = ldl (A, [ ]) ;

try
    [L, D, Parent, fl] = ldl (A, [1 2])
    ok = 0 ;
catch
    ok = 1 ;
end
if (~ok) error ('?'), end

try
    ldl
    ok = 0 ;
catch
    ok = 1 ;
end
if (~ok) error ('?'), end

try
    [L, D, Parent, fl] = ldl (1)
    ok = 0 ;
catch
    ok = 1 ;
end
if (~ok) error ('?'), end

try
    x = ldl (1,2) ;
    ok = 0 ;
catch
    ok = 1 ;
end
if (~ok) error ('?'), end

A =[ ...
1.7     0     0     0     0     0     0     0   .13     0
  0   1.0     0     0   .02     0     0     0     0   .01
  0     0   1.5     0     0     0     0     0     0     0
  0     0     0   1.1     0     0     0     0     0     0
  0   .02     0     0   2.6     0   .16   .09   .52   .53
  0     0     0     0     0   1.2     0     0     0     0
  0     0     0     0   .16     0   1.3     0     0   .56
  0     0     0     0   .09     0     0   1.6   .11     0
.13     0     0     0   .52     0     0   .11   1.4     0
  0   .01     0     0   .53     0   .56     0     0   3.1 ] ;
A = sparse (A) ;
b = [ ...
 .98 .64 .05
 .58 .20 .44
 .42 .37 .30
 .51 .78 .84
 .33 .68 .01
 .43 .46 .76
 .22 .56 .97
 .57 .79 .99
 .76 .05 .78
 .52 .60 .43 ] ;
P = [3 10 2 5 8 6 9 7 1 4] ;
I = speye (10) ;

[L, D, Parent, fl] = ldl (A) ;
err = norm ((L+I)*D*(L+I)'-A, 1) ;
fprintf ('err: %g fl: %g\n', err, fl) ;
if (err > 1e-14) error ('?'), end ;

Parent2 = etree (A) ;
if (any (Parent2 ~= Parent)) error ('?'), end

[L, D, Parent] = ldl (A) ;
err = norm ((L+I)*D*(L+I)'-A, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14) error ('?'), end ;

[L, D] = ldl (A) ;
err = norm ((L+I)*D*(L+I)'-A, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14) error ('?'), end ;

L2 = ldl (A) ;
err = norm (L - L2, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14) error ('?'), end ;

x = ldl (A, [ ], b) ;
err = norm (A*x-b, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14) error ('?'), end ;

[x, fl] = ldl (A, [ ], b) ;
err = norm (A*x-b, 1) ;
fprintf ('err: %g fl: %g\n', err, fl) ;
if (err > 1e-14) error ('?'), end ;

[L, D, Parent, fl] = ldl (A, P) ;
err = norm ((L+I)*D*(L+I)'-A(P,P), 1) ;
fprintf ('err: %g fl: %g\n', err, fl) ;
if (err > 1e-14) error ('?'), end ;

Parent2 = etree (A (P,P)) ;
if (any (Parent2 ~= Parent)) error ('?'), end

figure (1)
clf
subplot (2,2,1), spy (A),           title ('original matrix') ;
subplot (2,2,2), spy (A (P,P)),     title ('permuted matrix') ;
subplot (2,2,3), spy (L+D+L'),      title ('L+D+L''') ;
subplot (2,2,4), treeplot (Parent), title ('elimination tree') ;

[L, D, Parent] = ldl (A, P) ;
err = norm ((L+I)*D*(L+I)'-A(P,P), 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14) error ('?'), end ;

[L, D] = ldl (A, P) ;
err = norm ((L+I)*D*(L+I)'-A(P,P), 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14) error ('?'), end ;

L2 = ldl (A, P) ;
err = norm (L - L2, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14) error ('?'), end ;

x = ldl (A, P, b) ;
err = norm (A*x-b, 1) ;
fprintf ('err: %g\n', err) ;
if (err > 1e-14) error ('?'), end ;

[x, fl] = ldl (A, P, b) ;
err = norm (A*x-b, 1) ;
fprintf ('err: %g fl: %g\n', err, fl) ;
if (err > 1e-14) error ('?'), end ;

fprintf ('\nldl: all tests passed\n') ;
