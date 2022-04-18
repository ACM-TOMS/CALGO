% DEMO_VAR  call var1_likelihood.m with same parameters as demo_var.f90
%
%   DEMO_VAR FILENAME begins by reading values of a vector autoregressive
%   (VAR) time series, X, together with VAR(1) time series parameters A and
%   Sigma from the file FILENAME, which could be previously created with
%   INIT_DEMO_VAR. The function then uses var_ll from the vauto package to
%   calculate the function value of the corresponding VAR(1) likelihood
%   function and its gradient. After that it calls var1_likelihood to
%   obtain the same function value and gradients, asserts that the same
%   values are obtained and prints the result.

function demo_var(filename)
  addpath('../vauto/');
  if nargin==0, filename = '../demo/var_2_3.dat'; end
  [X, Sigma, A] = loadfile(filename);
  r = size(X,1);
  
  % CALCULATE THE LIKELIHOOD AND ITS GRADIENT WITH VAR_LL
  [ll, ok, lld] = var_ll(X, A, Sigma);
  assert(ok)
  Ai_vauto = reshape(lld(1:r*r), r, r);
  Sigi_vauto = ivech(lld(r*r+1:end));
  %
  [l, Aa, Siga] = var1_likelihood(X, A, Sigma);
  assert(almostequal(l, ll))
  assert(almostequal(Aa, Ai_vauto) && almostequal(Siga, Sigi_vauto));
  
  % PRINT LIKELIHOOD AND DERIVATIVE
  fprintf('l =\n  %.6f\n', ll);
  disp('Aa = ')
  for i=1:r, fprintf('  %9.6f', Aa(i,:)); fprintf('\n'); end
  disp('Siga = ')
  for i=1:r, fprintf('  %9.6f', Siga(i,:)); fprintf('\n'); end
  fprintf('\n');
end

function [X,Sig,A] = loadfile(fname)
  f = fopen(fname);
  r = fscanf(f,'%f',1);
  n = fscanf(f,'%f',1);
  X = fscanf(f,'%f', [n,r])';
  Sig = fscanf(f,'%f',[r,r])';
  A = fscanf(f,'%f',[r,r])';
  fclose(f);
end
