function [mu, A, lambda_0, kappa, theta_param, sigma, rho] ...
          = decode_pcsv_param(theta, decode)
%DECODE_PCSV_PARAM Decodes PCSV parameters from vector theta using decode.
%
%  [mu, A, lambda_0, kappa, theta_param, sigma, rho] ...
%    = decode_pcsv_param(theta, decode) decodes the PCSV model parameters
%    encoded in the vector theta. The function decode is expected to be
%    the function returned by encode_pcsv_param. Since a PCSV model may
%    contain NON-quadratic orthogonal matrices as parameters, the function
%    decode may have an individual state.
%
%    INPUT  theta: flat parameter vector to be decoded
%          decode: decoder function returned by encode_pcsv_param
%
%  See also ENCODE_PCSV_PARAM
%
% created by Benedikt Rudolph
% DATE: 16-Aug-2012

  p = decode(theta);
  mu = p(1).value;
  A = p(2).value;
  lambda_0 = p(3).value;
  kappa = p(4).value;
  theta_param = p(5).value;
  sigma = p(6).value;
  rho = p(7).value;
end
