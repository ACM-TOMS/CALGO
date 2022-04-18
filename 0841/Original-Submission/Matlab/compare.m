% This program compares results of the FORTRAN
%    and MATLAB versions of BHESS and some auxiliary
%    routines.  It assumes that the bhtest.f routine
%    and been run, generating mat.m
%    and that the matrix a, n and tol from mat.m
%    have been taken as input to abh3.m
%
%
%  Note that BackErr measures error internal to the 
%    matlab routine, where the rest of the functions
%    measure differences between matlab and Fortran executbles. 
%
comment =  'BackErr is usually 3 to 5 digits smaller ' 
comment = 'than the other quantities. This reflects the fact that'
comment = ' BackErr is an internally computed quantity.  '
comment = 'The other quantities reflect loss '
comment =' of accuracy due to the 14 digit accuracy of mat.m'
   z = inv(ltot); 
BackErr =  norm(a - z*b*inv(z))
   ForErr = norm(B-b) 
   NormZ1Err= norm(z1-z*ones(n,1))
   Norm1ZErr= norm(onez-ones(1,n)*z)
   Norm1InvZErr = norm(oneinvz-ones(1,n)*inv(z))
   CondZ = cond(z) 
comment= ' This is the 2-norm condition number'
