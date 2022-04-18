function [outlier, idx] = mahdist(x, varargin)
% MAHDIST Test for outliers in multivariate data using Mahalanobis distance
% and F-test 
%
% Syntax
%
% [outliers, idx] = MAHDIST(x, Name, Value)
%
% Input Arguments
%
% x                - Array of multivariate data
% (samples of each variate in separate row)
%
% Name-Value Pair Arguments
%
% Specify comma-separated pairs of Name, Value arguments. Name is the
% argument name and Value is the corresponding value. Name must appear
% inside single quotes (' '). You can specify up to 2 name and value
% pair arguments. Any order of the Name-Value pairs is allowed.
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
% MAHDIST(x1, x2, Name, Value) tests the input arrays x1 and x2 for
% outliers presence. Test supposes an approximately normal multivariate
% distribution of input data.
%
% Example
%
% Test the input bivariate array inputArray for outliers presence with
% verbose output.
%
% inputArray = [154 136 91 125 133 125 93 80 132 107 142 115 114 120 141; 108  90 54  89  93  77 43 50 125  76  96  74  79  71  90];
% [outlier, idx] = mahdist(inputArray, 'verboseOutput', 'on');
%
% Outliers are pair 9 (132 125) and pair 7 (93 43).
%
% More About
%
% For reference see ch. 5.1 of
% Afifi, A.A., Azen, S.P. 1979. Statistical analysis: a computer oriented
% approach (2nd ed.). Academic Press.
% Example data are taken from this source.
%
% See also MAHAL, FCDF.
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

for curArg = 1 : 2 : nargin-2
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

xCopy = x;

outlier = [];
idx = [];

if verboseOutput
    fprintf('Significance level: %.2f\n', alpha);
    fprintf('Variates:           %d\n', size(x, 1));
    fprintf('Samples:            %d\n', size(x, 2));
end

while size(x, 2) > 2
    [minP, idx_] = findOutlier(xCopy, verboseOutput);

    if minP < alpha
        outlier_ = xCopy(:, idx_);
        outlier = horzcat(outlier, outlier_);

        [~, idxResult] = ismember(xCopy(:, idx_)', x', 'rows');
        idx = [idx idxResult];

        xCopy(:, idx_) = [];
    else
        break
    end
end

%% find the outlier
% input:  x1   - vector of variate 1
%         x2   - vector of variate 2
%         verboseOutput - verbose output (true or false)
% output: minP - min P
%         idx  - index of potential outlier
function [minP, idx] = findOutlier(x, verboseOutput)

n = size(x, 2);
k = n - 1;
p = size(x, 1);

P = zeros(1, n);

if verboseOutput
    fprintf('----------------------------------------\n');
    fprintf('Sample Test statistic P\n');
end

for i = 1 : n
    cand = x(:, i);

    x_ = x;
    x_(:, i) = [];

    D = mahal(cand', x_');

    F = (k - p) * k * D / (k^2 - 1) / p;
    P(i) = 1 - fcdf(F, p, k-p);

    if verboseOutput
        fprintf('%4d\t%6.2f\t\t%f\n', i, F, P(i));
    end
end

[minP, idx] = min(P);

if verboseOutput
    fprintf('min P: %f, Sample: %d (', minP, idx);

    for i = 1 : size(x, 1)
        fprintf('%.1f ', x(i, idx));
    end

    fprintf(')\n');
end

end

end