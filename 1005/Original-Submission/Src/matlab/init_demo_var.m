% INIT_DEMO_VAR  Initialize data files for demo_var
%
%   INIT_DEMO_VAR(R,N) generates a random vector-autoregressive VAR(1) time
%   series of dimension R and length N, and writes it to a file named
%   var_<R>_<N>.dat, for example var_3_50.dat for R=3 and N=50. The covariance
%   matrix Sigma and the coefficent matrix A are also written to this file.
%
%   The program varma_sim from [1] is used to generate the series, using the
%   r x r Lehmer matrix for Sigma and a random matrix for A.
%
%   [1] Jonasson K: Algorithm 878: Exact VARMA likelihood and its gradient for complete and
%       incomplete data with Matlab, ACM TOMS 35, 2008.

function init_demo_var(r, n)
  addpath('../vauto');
  rand('twister', 123645) %#ok<RAND>
  Sigma = gallery('lehmer', r);
  A = rand(r,r) .* 1/(r*r);
  X = varma_sim(A, [], Sigma, n); %#ok<NASGU>
  name = ['var_' int2str(r) '_' int2str(n) '.dat'];
  f = fopen(name, 'w');
  fprintf(f, '%d\n%d\n', r, n);
  fclose(f);
  save(name, 'X', 'Sigma', 'A', '-ascii', '-append', '-double')
end
