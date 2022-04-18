function [outliers, idx] = gesd(x, varargin)
% GESD Two-sided Generalized (extreme Studentized deviate) ESD test for one
% or more outliers
%
% Syntax
%
% [outliers, idx] = GESD(x, Name, Value)
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
% 'outliersNumber' - Upper bound for suspected number of outliers (mandatory)
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
% GESD(x, Name, Value) tests the input array x for one or more outliers
% presence. Test supposes an approximately normal distribution of input data.
%
% Example
%
% Test the input array inputArray for presence of 10 outliers with verbose
% output.
%
% inputArray = [-0.25 0.68 0.94 1.15 1.20 1.26 1.26 1.34 1.38 1.43 1.49 1.49 1.55 1.56 1.58 1.65 1.69 1.70 1.76 1.77 1.81 1.91 1.94 1.96 1.99 2.06 2.09 2.10 2.14 2.15 2.23 2.24 2.26 2.35 2.37 2.40 2.47 2.54 2.62 2.64 2.90 2.92 2.92 2.93 3.21 3.26 3.30 3.59 3.68 4.30 4.64 5.34 5.42 6.01];
% [outliers, idx] = gesd(inputArray, 'outliersNumber', 10, 'verboseOutput', 'on');
%
% More About
%
% For reference see NIST/SEMATECH e-Handbook of Statistical Methods,
% http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm
% (example data are taken from this source).
%
% See also DIXON, GRUBBS, MZSCORE, TIETJEN, TUKEY.
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

R        = zeros(1, k);
Lambda   = zeros(1, k);
outliers = zeros(1, k);

n = length(x);
xCopy = x;

for i = 1 : k
    xMean = mean(x);
    xSD = std(x);

    R(i) = max( abs(x - xMean) / xSD );

    % compute critical value
    p = 1 - alpha /2 /(n - i + 1);
    t = tinv(p, n - i - 1);
    Lambda(i) = t * (n - i) / sqrt((n - i - 1 + t^2) * (n - i + 1));

    if R(i) > Lambda(i)
        outliersNumber = i;
    end

    [~, idxMax] = max( abs(x - xMean) );
    outliers(i) = x(idxMax);
    x(idxMax) = [];
end

idx = zeros(1, outliersNumber);

for i = 1 : outliersNumber
    idx(i) = find(xCopy == outliers(i), 1);
    xCopy( idx(i) ) = NaN;
end

if outliersNumber > 0
    outliers = outliers(1 : outliersNumber);
end

if verboseOutput
    fprintf('Significance level: %.2f\n', alpha);

    for i = 1 : k
        fprintf('Outliers: %2d, Test statistic: %.3f, critical value: %.3f\n', ...
            i, R(i), Lambda(i));
    end

    fprintf('Number of outliers: %d\n', outliersNumber);
end

end