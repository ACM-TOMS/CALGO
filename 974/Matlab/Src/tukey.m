function [outliers, idx] = tukey(x, varargin)
% TUKEY Tukey test to screen data for outliers
%
% Syntax
%
% [outliers, idx] = TUKEY(x, Name, Value)
%
% Input Arguments
%
% x                - Input array
%
% Name-Value Pair Arguments
%
% Specify comma-separated pairs of Name, Value arguments. Name is the
% argument name and Value is the corresponding value. Name must appear
% inside single quotes (' '). You can specify up to 2 name and value
% pair arguments. Any order of the Name-Value pairs is allowed.
%
% 'whisker'        — Maximum whisker length (optional)
% 1.5 (default) | positive numeric value
%
% 'plotBoxplot'    - Boxplot is plotted (optional)
% true | false (default)
%
% Output Arguments
%
% outliers         - Outliers
% idx              - Indices of outliers in input array
%
% Description
%
% TUKEY(x, Name, Value) screens the input array x for multiple outliers
% presence.
%
% Example
%
% Test the input array inputArray for outlier presence and plot boxplot.
%
% [outliers, idx] = tukey(inputArray, 'plotBoxplot', true);
%
% See also BOXPLOT, DIXON, GRUBBS, GESD, MZSCORE, TIETJEN.
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
whisker = 1.5;
plotBoxplot = false;

okParName = {'whisker', 'plotBoxplot'};

for curArg = 1 : 2 : nargin-1
    parName  = varargin{curArg};
    parValue = varargin{curArg + 1};

    parID = find(strcmpi(parName, okParName));

    if isempty(parID)
        error( ['Wrong input argument: Unknown parameter ''', parName, ''''] )
    else
        switch(parID)
            case 1 % whisker

                whisker = parValue;

            case 2 % plotBoxplot

                switch parValue
                    case {true, false}
                        plotBoxplot = parValue;
                    otherwise
                        error('Wrong input argument: plotBoxplot should be true or false');
                end
        end
    end
end

%% test calculation

Q1 = quantile(x, 0.25);
Q3 = quantile(x, 0.75);
IQR = Q3 - Q1;  % interquartile range

lowest  = Q1 - whisker * IQR;
highest = Q3 + whisker * IQR;

idxLow  = find(x < lowest);
idxHigh = find(x > highest);

idx = [idxLow idxHigh];
outliers = x(idx);

if plotBoxplot
    boxplot(x, 'whisker', whisker, 'symbol', 'r.')
end

end