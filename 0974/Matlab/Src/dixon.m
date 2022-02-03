function [outlier, idx] = dixon(x, varargin)
% DIXON One-sided Dixon test for single outlier
%
% Syntax
%
% [outliers, idx] = DIXON(x, Name, Value)
%
% Input Arguments
%
% x             - Input array (size should be between 3 and 30)
%
% Name-Value Pair Arguments
%
% Specify comma-separated pairs of Name, Value arguments. Name is the
% argument name and Value is the corresponding value. Name must appear
% inside single quotes (' '). You can specify up to 2 name and value
% pair arguments. Any order of the Name-Value pairs is allowed.
%
% 'alpha'       - Significance level (optional)
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
% DIXON(x, Name, Value) tests the input array x for single outlier
% presence. Test supposes an approximately normal distribution of input data.
%
% Example
%
% Test the input array inputArray for outlier presence with
% significance level 0.10 and verbose output.
%
% inputArray = [568 570 570 570 572 578 584 596];
% [outlier, idx] = dixon(inputArray, 'verboseOutput', 'on', 'alpha', 0.1);
%
% More About
%
% For reference see the documentation of Dataplot software developed by NIST,
% http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/dixon.htm
% (example data are taken from this source).
%
% See also GESD, GRUBBS, MZSCORE, TIETJEN, TUKEY.
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

%% constants

SIDE_LEFT  = 'left';
SIDE_RIGHT = 'right';

simulations = 25000;

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

n = length(x);

if (n < 3) || (n > 30)
    error('Wrong input argument: size of input array should be between 3 and 30');
end

%% test calculation

[D, side] = getTestStatistic(x, n);

Dcrit = getCriticalValue(n, simulations, alpha);

if D > Dcrit
    xSort_ = sort(x);

    if strcmp(side, SIDE_LEFT)
        outlier = xSort_(1);
    else
        outlier = xSort_(end);
    end

    idx = find(x == outlier, 1);
else
    outlier = [];
    idx = [];
end

if verboseOutput
    fprintf('Significance level: %.2f\n',    alpha);
    fprintf('Test kind:          %s-side\n', side);
    fprintf('Test statistic:     %f\n',      D);
    fprintf('Critical value:     %f (based on %d simulations)\n', Dcrit, simulations);
end

%% compute critical value by simulation
% input:  n           - number of observations
%         simulations - number of simulations
%         alpha       - significance level
% output: Dcrit       - critical value
function Dcrit = getCriticalValue(n, simulations, alpha)

dSimul = zeros(simulations, 1);
pd = makedist('Normal');

case_3_7   = num2cell( 3 :  7);
case_8_10  = num2cell( 8 : 10);
case_11_13 = num2cell(11 : 13);
case_14_30 = num2cell(14 : 30);

for i = 1 : simulations
    x_ = random(pd, 1, n);
    xSort = sort(x_);

    switch(n)
        case case_3_7      % test r10
            dSimul(i) = (xSort(2) - xSort(1)) / (xSort(n) - xSort(1));

        case case_8_10     % test r11
            dSimul(i) = (xSort(2) - xSort(1)) / (xSort(n-1) - xSort(1));

        case case_11_13    % test r21
            dSimul(i) = (xSort(3) - xSort(1)) / (xSort(n-1) - xSort(1));

        case case_14_30    % test r22
            dSimul(i) = (xSort(3) - xSort(1)) / (xSort(n-2) - xSort(1));
    end
end

Dcrit = quantile(dSimul, 1 - alpha);

end

%% get the test statistic
% input:  x    - vector of observations
%         n    - number of observations
% output: D    - test statistic
%         side - test kind (left-side or right-side)
function [D, side] = getTestStatistic(x, n)

% find more extreme between minimum and maximum points

xSort = sort(x);
difMin = xSort(2) - xSort(1);
difMax = xSort(end) - xSort(end-1);

if difMin > difMax
    side = SIDE_LEFT;
else
    side = SIDE_RIGHT;
end

case_3_7   = num2cell( 3 :  7);
case_8_10  = num2cell( 8 : 10);
case_11_13 = num2cell(11 : 13);
case_14_30 = num2cell(14 : 30);

if strcmp(side, SIDE_LEFT) % test for minimum
    switch(n)
        case case_3_7      % test r10
            D = (xSort(2) - xSort(1)) / (xSort(n) - xSort(1));

        case case_8_10     % test r11
            D = (xSort(2) - xSort(1)) / (xSort(n-1) - xSort(1));

        case case_11_13    % test r21
            D = (xSort(3) - xSort(1)) / (xSort(n-1) - xSort(1));

        case case_14_30    % test r22
            D = (xSort(3) - xSort(1)) / (xSort(n-2) - xSort(1));
    end
else                       % test for maximum
    switch(n)
        case case_3_7      % test r10
            D = (xSort(n) - xSort(n-1)) / (xSort(n) - xSort(1));

        case case_8_10     % test r11
            D = (xSort(n) - xSort(n-1)) / (xSort(n) - xSort(2));

        case case_11_13    % test r21
            D = (xSort(n) - xSort(n-2)) / (xSort(n) - xSort(2));

        case case_14_30    % test r22
            D = (xSort(n) - xSort(n-2)) / (xSort(n) - xSort(3));
    end
end

end

end