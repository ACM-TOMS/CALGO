function [outliers, idx] = tietjen(x, varargin)
% TIETJEN Two-sided Tietjen-Moore test for multiple outliers
%
% Syntax
%
% [outliers, idx] = TIETJEN(x, Name, Value)
%
% Input Arguments
%
% x                - Input array
%
% Name-Value Pair Arguments
%
% Specify comma-separated pairs of Name, Value arguments. Name is the
% argument name and Value is the corresponding value. Name must appear
% inside single quotes (' '). You can specify up to 3 name and value
% pair arguments. Any order of the Name-Value pairs is allowed.
%
% 'outliersNumber' - Number of outliers (mandatory)
%
% 'alpha'          - Significance level (optional)
% value between 0 and 1
% if not provided, default is 0.05 for 5% significance
%
% 'verboseOutput'  - Verbose output (optional)
% 'on' | 'off' (default)
%
% Output Arguments
%
% outliers         - Outliers
% idx              - Indices of outliers in input array
%
% Description
%
% TIETJEN(x, Name, Value) tests the input array x for multiple outliers
% presence. Test supposes an approximately normal distribution of input data.
%
% Example
%
% Test the input array inputArray for presence of 2 outliers with verbose
% output.
%
% inputArray = [-1.40 -0.44 -0.30 -0.24 -0.22 -0.13 -0.05 0.06 0.10 0.18 0.20 0.39 0.48 0.63 1.01];
% [outliers, idx] = tietjen(inputArray, 'outliersNumber', 2, 'verboseOutput', 'on');
%
% More About
%
% For reference see NIST/SEMATECH e-Handbook of Statistical Methods,
% http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h2.htm
% (example data are taken from this source).
%
% See also DIXON, GRUBBS, GESD, MZSCORE, TUKEY.
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

% number of input arguments must be 3-7
narginchk(3, 7)

% default values of input arguments
alpha = 0.05;
verboseOutput = false;

okParName = {'alpha', 'verboseOutput', 'outliersNumber'};

for curArg = 1 : 2 : nargin-1
    parName  = varargin{curArg};
    parValue = varargin{curArg + 1};

    parID = find(strcmpi(parName, okParName));

    if isempty(parID)
        error( ['Wrong input argument: Unknown parameter ''', parName, ''''] )
    else
        switch(parID)
            case 1 % alpha
                alpha = parValue;

            case 2 % verboseOutput

                switch parValue
                    case 'on'
                        verboseOutput = true;
                    case 'off'
                        verboseOutput = false;
                    otherwise
                        error('Wrong input argument: verbose output should be ''on'' or ''off'' ');
                end

            case 3 % outliersNumber
                k = parValue;
        end
    end
end

%% test calculation

n = length(x);

[Ek, idx] = getTestStatistic(x, n, k);

% compute critical value by simulation

simulations = 10000;
ekSimul = zeros(simulations, 1);
pd = makedist('Normal');

for i = 1 : simulations
    x_ = random(pd, 1, n);
    [ekSimul(i), ~] = getTestStatistic(x_, n, k);
end

Ecrit = quantile(ekSimul, alpha);

if Ek < Ecrit
    idx = sort(idx);
    outliers = x(idx);
else
    outliers = [];
    idx = [];
end

if verboseOutput
    fprintf('Significance level: %.2f\n', alpha);
    fprintf('Test statistic:     %f\n', Ek);
    fprintf('Critical value:     %f (based on %d simulations)\n', Ecrit, simulations);
end

%% get the test statistic
% input:  x   - vector of observations
%         n   - number of observations
%         k   - number of outliers
% output: Ek  - test statistic
%         idx - indices of k outliers
function [Ek, idx] = getTestStatistic(x, n, k)

% compute the absolute residuals
r = abs(x - mean(x));

% sort the input values by their absolute residuals in ascending order
[~, resOrder] = sort(r);
z = x(resOrder);

% get the indices of k candidates for outliers
idx = resOrder(n-k+1 : end);

% create the data subset without k largest values
subsetX = z(1 : n-k);

% compute the sums of squares
sumSubset = (subsetX - mean(subsetX)) .^ 2;
sumAll = (x - mean(x)) .^ 2;

% compute the test statistic
Ek = sum(sumSubset) / sum(sumAll);

end

end