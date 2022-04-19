function [outliers, idx] = mzscore(x)
% MZSCORE Modified Z-score to screen data for outliers 
%
% Syntax
%
% [outliers, idx] = MZSCORE(x)
%
% Input Arguments
%
% x        - Input array
%
% Output Arguments
%
% outliers - Potential outliers
% idx      - Indices of potential outliers in input array
%
% Description
%
% MZSCORE(x) screens the input array x for multiple outliers presence.
% Test supposes an approximately normal distribution of input data.
%
% More About
%
% For reference see NIST/SEMATECH e-Handbook of Statistical Methods,
% http://itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
%
% See also DIXON, GESD, GRUBBS, TIETJEN, TUKEY.
%
% Information about MATLAB version, revision and date
%
% MATLAB version: R2015a
% Revision:       1.0
% Date:           2016/05/04
%
% Authors
%
% A. Novoselsky, The Weizmann Institute of Science
% E. Kagan, Ariel University

%% check of input arguments

if nargin == 0
    error('Input arguments are NOT passed to function')
end

% number of input arguments must be 1
narginchk(1, 1)

%% test calculation

n = length(x);
M = zeros(1, n);

medianX = median(x);
MAD = median( abs(x - medianX) );

for i = 1 : n
    M(i) = 0.6745 * abs(x(i) - medianX) / MAD;
end

idx = find(M > 3.5);
outliers = x(idx);

end