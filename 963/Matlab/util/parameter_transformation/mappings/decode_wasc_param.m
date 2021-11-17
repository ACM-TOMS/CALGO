function [mu, Sigma_0, M, Q, rho, beta] = decode_wasc_param(theta, decode)
%DECODE_WASC_PARAM Decodes WASC parameters from vector theta using decode.
%
%  [mu, Sigma_0, M, Q, rho, beta] = decode_wasc_param(theta, decode)
%    decodes the PCSV model parameters encoded in the vector theta.
%    The function decode is expected to be the function returned by
%    encode_wasc_param.
%
%    INPUT  theta: flat parameter vector to be decoded
%          decode: decoder function returned by encode_wasc_param
%
%  See also ENCODE_WASC_PARAM
%
% created by Benedikt Rudolph
% DATE: 16-Aug-2012

  p = decode(theta);
  mu = p(1).value;
  Sigma_0 = p(2).value;
  M = p(3).value;
  Q = p(4).value;
  rho = p(5).value;
  beta = p(6).value;
end
