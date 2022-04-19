function [s, exitflag, fNotEval, varNotUsed] = compSigma(fcn, n, varargin)


global DAEZERO; DAEZERO = sigma(n);
exitflag = 0; % initialize exitflag
s = zeros(n); % initialize sigma

y(1:n) = sigma(n);
t = sigma(n);
for i = 1:n
    y(i) = set(y(i), i);
end

f = feval(fcn, t, y, varargin{:});
% if the number of variables is larger than problem size, it's a Matlab error

m = length(f); % number of equations
if m>n % more equations than needed
    s = [];
    exitflag = -3; % cannot visualize sigma, at least for now
    return
end

fNotEval = [];
for i = 1:m
    if isa(f(i), 'numeric')
        s(i,:) = repmat(-inf, [1 n]);
        exitflag = -2; % ill-posed, but still can visualize sigma
        fNotEval = [fNotEval i];
    else
        s(i,:) = getVector(f(i));
        if any(isfinite(s(i,:)))==0
            exitflag = -2; % ill-posed, but still can visualize sigma
            fNotEval = [fNotEval i];
        end
    end
end

% Less equations
if m<n
    row = (m+1):n;
    s(row,:) = repmat(-inf, [n-m, n]);
    exitflag = -2; % ill-posed, but still can visualize sigma
    fNotEval = [fNotEval row];
end

if isempty(fNotEval)
    fNotEval = [];
end

% Check variables
checkVar = any(isfinite(s), 1);
varNotUsed = find(checkVar==0);
if ~isempty(varNotUsed)
    exitflag = -2; % ill-posed, but still can visualize sigma
else
    varNotUsed = [];
end
end