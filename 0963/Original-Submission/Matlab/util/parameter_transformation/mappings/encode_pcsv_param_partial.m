function [theta, decode] = encode_pcsv_param_partial(mu, A, lambda_0, kappa ...
                                              , theta, sigma, rho)
%ENCODE_PCSV_PARAM_PARTIAL Encodes all parameters of a PCSV partial model
%
%  [theta, decode] = ENCODE_PCSV_PARAM_PARTIAL(mu, A, lambda_0, kappa, theta, sigma, 
%    rho) encodes the parameters mu, A, lambda_0, kappa, theta, sigma, rho of
%    a PCSV partial model into a single parameter vector theta. The returned function
%    decode can be used together with decode_pcsv_param to recover the
%    parameters from theta.
%
% See also DECODE_PCSV_PARAM, ENCODE_PARAMETERS, ENCODE_WASC_PARAM.
%
% created by Benedikt Rudolph
% DATE: 30-Jan-2013

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
  param(end).value = kappa;

  param(end+1).name = 'theta';
  param(end).constraint = 'none';
  param(end).value = theta;

  param(end+1).name = 'sigma';
  param(end).constraint = 'none';
  %param(end).constraint = 'equal';
  param(end).value = sigma;

  param(end+1).name = 'rho';
  param(end).constraint = 'none';
  param(end).value = rho;
  
  [theta, decode] = encode_parameters(param);
end
