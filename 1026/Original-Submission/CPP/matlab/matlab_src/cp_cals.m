function [U,varargout] = cp_cals(X,R,ktensors,varargin)

%  cp_cals Compute multiple CP decompositions using CALS
%
%  It works similarly to, and is mostly compatible with, the CP_ALS
%  function provided by the Tensor Toolbox. See CP_ALS for more information.
%
%  The differences with CP_ALS are:
%    - The input tensor must be a dense tensor
%    - The 'dimorder' parameter is not accepted
%
%  Usage:
%
%  M = cp_cals(X,R,ktensors) computes the CP decomposition of len(R) models using an alternating least-squares algorithm.
%     'X' - The target tensor.
%	  'R' - Vector of ranks of the ktensors to fit to the target tensor.
%	  'ktensors' - Either 'random', to use randomly initialized ktensors, or a cell array of len(R) ktensors.
%	The result M is a cell array of ktensors.
%
%  [M,U0] = cp_cals(...) also returns the initial guesses.
%
%  M = cp_cals(X,R,ktensors,'param',value,...) specifies optional parameters and values.
%
%      Valid parameters and their default values are:
%        'tol' - Tolerance on difference in fit {1.0e-4}
%        'maxiters' - Maximum number of iterations {50}
%        'buffer-size' - Maximum size of the factor matrices buffer {4200}
%        'update-method' - Method for updating the factor matrices {unconstrained}
%                          ['unconstrained', 'nnls']
%        'mttkrp-method' - method for computing MTTKRP
%                          ['mttkrp', 'twostep0', 'twostep1', 'auto']
%
%        'no-ls'(default)/'ls' - Whether to use line search.
%        'ls-interval' - Interval (per individual model) to apply line search {5}
%        'ls-step' - Factor with which to jump in line search {1.2}
%
%        'no-cuda'(default)/'cuda' - Whether to use cuda (make sure the binary is compiled with cuda support)

argout = cell(1, nargout-1);
[U,argout{:}] = cp_cals_driver(X,R,ktensors,varargin{:});
varargout = argout;
