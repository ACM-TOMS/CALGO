function [theta, decode] = encode_wasc_param(mu, Sigma_0, M, Q, rho, beta)
%ENCODE_WASC_PARAM Encodes all parameters of a WASC model into a single vector
%
%  [theta, decode] = ENCODE_WASC_PARAM(mu, Sigma_0, M, Q, rho, beta) encodes
%    the parameters mu, Sigma_0, M, Q, rho, beta of a WASC model into a single
%    parameter vector theta. The returned function decode can be used together
%    with decode_wasc_param to recover the parameters from theta. Since a WASC
%    model contains no NON-quadratic orthogonal matrices as parameters, the
%    returned function decode does not have an individual state. Hence, any
%    returned decode function from encode_wasc_param can be used to recover
%    the parameters (as opposed to PCSV models).
%
% See also DECODE_WASC_PARAM, ENCODE_PARAMETERS, ENCODE_PCSV_PARAM.
%
% created by Benedikt Rudolph
% DATE: 16-Aug-2012

  param(1).name = 'mu';
  param(end).constraint = 'equal';
  param(end).value = mu;

  param(end+1).name = 'Sigma_0';
  %param(end).constraint = 'spd';
  param(end).constraint = 'equal';
  param(end).value = Sigma_0;

  param(end+1).name = 'M';
  param(end).constraint = 'None';
  param(end).value = M;

  param(end+1).name = 'Q';
  param(end).constraint = 'None';
  param(end).value = Q;

  param(end+1).name = 'rho';
  param(end).constraint = 'None';
  param(end).value = rho;

  param(end+1).name = 'beta';
  param(end).constraint = 'None';
  param(end).value = beta;
  
  [theta, decode] = encode_parameters(param);
end
