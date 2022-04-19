%*****************************************************************************
% DSDP5:  Dual-Scaling Algorithm for Positive Semidefinite Programming
% Copyright (c) 2004 by
% S. J. Benson, Y. Ye
% Last modified: 20 October 2004
%*****************************************************************************
% check;
%
% This script tests the DSDP solver using the examples. Compare the output
% to the file check.out
% 
diary on;
ntrials=2;
BB=cell(1);
for trials = [1:ntrials];
  N = 10;
  B = zeros(N); 
  for i=1:N,
    for j=1:N, B(i,j)=mod(i+j,trials+1); end;
  end;
  fprintf('MAXCUT \n');
  [y,X,obj] = maxcut(B); 
  fprintf('GPP \n');
  [y,X,obj] = gpp(B,1); 
  fprintf('ETP \n');
  [obj,d]=etp(B*B'+speye(N)');
  fprintf('THETA\n');
  [obj,X]=thetaproblem(B);
end;

%[X,y] = randdinf([10 4 3], 20);

OPTIONS=doptions;
OPTIONS.print=10;

fprintf('CONTROL1\n');
[AC,b] = readsdpa('control1.dat-s');
[STAT,y,X]=dsdp(b,AC,OPTIONS);

fprintf('ARCH0\n');
[AC,b] = readsdpa('arch0.dat-s');
[STAT,y,X]=dsdp(b,AC,OPTIONS);

diary off;
