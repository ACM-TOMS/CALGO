function [theta, decode] = encode_parameters(param)
%ENCODE_PARAMETERS Encodes the struct array param into a single vector

%  [theta, decode] = encode_parameters(param) encodes the list of parameters
%    in param into a flat parameter vector theta  which can be decoded by the
%    returned function decode.
%    param is expected to be a struct array with the following fields:
%    - param.name the parameter's name
%    - param.constraint a string representing the constraint to be
%      enforced on this parameter. Can be one of
%      'equal', 'positive', 'correlation', 'orthogonal_quadratic',
%      'orthogonal_non_quadratic', 'spd'.
%       Any other value will be considered as no constraint.
%    - param.value the value to be encoded
%    For individual models there exist parameter encoding and decoding
%    functions that conveniently wrap away the setup of param.
%
% See also ENCODE_PCSV_PARAM, ENCODE_WASC_PARAM.
%
% created by Benedikt Rudolph
% DATE: 20-Aug-2012
  
  theta = [];
  
  for j=1:length(param) % iterate over all parameters
    v = param(j).value;
    param(j).size = size(v);
    last_idx = length(theta);
    switch param(j).constraint
      case 'equal' % forced to equal initial value
        param(j).index = [];
      case 'positive'
        v = pos_to_real(v(:))';
        theta = [theta, v];
        param(j).index = last_idx + (1:length(v));
      case 'correlation'
        v = corr_to_real(v(:))';
        theta = [theta, v];
        param(j).index = last_idx + (1:length(v));
      case 'orthogonal_quadratic'
        % v is assumed to be quadratic
        p = size(v,2);
        % perform a Cayley transformation of Q to obtain
        % the corresponding skew-symmetric matrix S
        S = (eye(p) - v) * inv(eye(p) + v);
        s = S(logical(tril(ones(p),-1)))';
        theta = [theta, s];
        param(j).index = last_idx + (1:length(s));
      case 'orthogonal_non_quadratic'
        % v can be non-quadratic
        p = size(v,2);
        % actual base value is stored in param(j).value
        % theta encodes skew symmetric matrix for Cayley transformation
        % initially that is just zeros
        % see decode function for details
        theta = [theta, zeros(1,p*(p-1)/2)];
        param(j).index = last_idx + (1:p*(p-1)/2);
      case 'spd' % symmetric positive definite
        p = size(v,1);
        R = chol(v);
        r = R(logical(triu(ones(p))))';
        theta = [theta, r];
        param(j).index = last_idx + (1:length(r));
      otherwise % no constraint
        v = v(:)';
        theta = [theta, v];
        param(j).index = last_idx + (1:length(v));
    end
  end
  
  decode = @(th) decode_parameters(th, param);
end

function param = decode_parameters(theta, p0)
  for j=1:length(p0) % iterate over all parameters
    param(j).name = p0(j).name; % copy name
    th = theta(p0(j).index); % extract relevant part of theta
    sz = p0(j).size;
    switch p0(j).constraint
      case 'equal' % forced to equal initial value
        param(j).value = p0(j).value;
      case 'positive'
        param(j).value = reshape(real_to_pos(th), sz);
      case 'correlation'
        param(j).value = reshape(real_to_corr(th), sz);
      case 'orthogonal_quadratic'
        p = sz(2);
        % obtain skew symmetric matrix
        S = vector_to_skew_symmetric(th, p);
        % perform Cayley transformation
        param(j).value = (eye(p) + S) \ (eye(p) - S);
      case 'orthogonal_non_quadratic'
        A0 = p0(j).value;
        p = sz(2);
        % obtain skew symmetric matrix
        S = vector_to_skew_symmetric(th, p);
        % perform Cayley transformation
        Q = (eye(p) + S) \ (eye(p) - S);
        param(j).value = A0 * Q;
      case 'spd'
        p = sz(1);
        R = zeros(p);
        R(logical(triu(ones(p)))) = th;
        param(j).value = R'*R;
      otherwise % no constraint
        param(j).value = reshape(th, sz);
    end
  end
end
