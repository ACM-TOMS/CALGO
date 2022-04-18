function [theta, decode] = encode_pcsv_param(mu, A, lambda_0, kappa ...
                                              , theta, sigma, rho)
%ENCODE_PCSV_PARAM Encodes all parameters of a PCSV model into a single vector
%
%  [theta, decode] = ENCODE_PCSV_PARAM(mu, A, lambda_0, kappa, theta, sigma, 
%    rho) encodes the parameters mu, A, lambda_0, kappa, theta, sigma, rho of
%    a WASC model into a single parameter vector theta. The returned function
%    decode can be used together with decode_pcsv_param to recover the
%    parameters from theta. A PCSV model may contain NON-quadratic orthogonal
%    matrices as parameter and hence the returned function decode has and
%    individual state. That means that two returned decode function from two
%    calls to encode_pcsv_param are different and cannot be interchanged.
%    If and only if the numbers of eigenvectors is equal to the number of
%    assets (n==p), the returned function decode has no individual state.
%
% See also DECODE_PCSV_PARAM, ENCODE_PARAMETERS, ENCODE_WASC_PARAM.
%
% created by Benedikt Rudolph
% DATE: 16-Aug-2012

  [n, p] = size(A);
  
  param(1).name = 'mu';
  %param(end).constraint = 'None';
  param(end).constraint = 'equal';
  param(end).value = reshape(mu, n, 1);

  param(end+1).name = 'A';
  %param(end).constraint = 'orthogonal_non_quadratic';
  %% use more flexible constraint trafo for quadratic orthogonal matrices
  %if n==p && det(A)>0
    %param(end).constraint = 'orthogonal_quadratic';
  %end
  param(end).constraint = 'equal';
  param(end).value = A;

  param(end+1).name = 'lambda_0';
  param(end).constraint = 'equal';
  param(end).value = reshape(lambda_0, 1, p);

  param(end+1).name = 'kappa';
  param(end).constraint = 'none';
  param(end).value = reshape(kappa, 1, p);

  param(end+1).name = 'theta';
  param(end).constraint = 'none';
  param(end).value = reshape(theta, 1, p);

  param(end+1).name = 'sigma';
  param(end).constraint = 'none';
  %param(end).constraint = 'equal';
  param(end).value = reshape(sigma, 1, p);

  param(end+1).name = 'rho';
  param(end).constraint = 'none';
  param(end).value = reshape(rho, 1, p);
  
  [theta, decode] = encode_parameters(param);
end
