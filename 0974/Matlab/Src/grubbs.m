function [outlier, idx] = grubbs(x, varargin)
% GRUBBS Two-sided Grubbs test for single outlier
%
% Syntax
%
% [outliers, idx] = GRUBBS(x, Name, Value)
%
% Input Arguments
%
% x             - Input array
%
% Name-Value Pair Arguments
%
% Specify comma-separated pairs of Name, Value arguments. Name is the
% argument name and Value is the corresponding value. Name must appear
% inside single quotes (' '). You can specify up to 2 name and value
% pair arguments. Any order of the Name-Value pairs is allowed.
%
% 'alpha'        - Significance level (optional)
% value between 0 and 1
% if not provided, default is 0.05 for 5% significance
%
% 'verboseOutput' - Verbose output (optional)
% 'on' | 'off' (default)
%
% Output Arguments
%
% outlier       - Outlier
% idx           - Index of outlier in input array
%
% Description
%
% GRUBBS(x, Name, Value) tests the input array x for single outlier
% presence. Test supposes an approximately normal distribution of input data.
%
% Example
%
% Test the input array inputArray for outlier presence with
% significance level 0.10 and verbose output.
%
% inputArray = [199.31, 199.53, 200.19, 200.82, 201.92, 201.95, 202.18, 245.57];
% [outlier, idx] = grubbs(inputArray, 'verboseOutput', 'on', 'alpha', 0.1);
%
% More About
%
% Grubbs test is also known as the maximum normed residual test.
% For reference see NIST/SEMATECH e-Handbook of Statistical Methods,
% http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h1.htm
% (example data are taken from this source).
%
% See also DIXON, GESD, MZSCORE, TIETJEN, TUKEY.
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

% number of input arguments must be 1-5
narginchk(1, 5)

% default values of input arguments
alpha = 0.05;
verboseOutput = false;

okParName = {'alpha', 'verboseOutput'};

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
        end
    end
end

%% test calculation

n = length(x);
xMean = mean(x);
xSD = std(x);

% Grubbs test statistic
G = max( abs(x - xMean) / xSD );

tCritValue = tinv(alpha /2 /n, n - 2);
Gcrit = (n - 1) / sqrt(n) * sqrt( tCritValue^2 / (n - 2 + tCritValue^2) );

if G > Gcrit
    [~, idx] = max( abs(x - xMean) );
    outlier = x(idx);
else
    outlier = [];
    idx = [];
end

if verboseOutput
    fprintf('Significance level: %.2f\n', alpha);
    fprintf('Test statistic:     %f\n', G);
    fprintf('Critical value:     %f\n', Gcrit);
end

end