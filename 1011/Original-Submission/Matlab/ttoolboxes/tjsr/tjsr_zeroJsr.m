function [out,normHist] = tjsr_zeroJsr(M,tol)
%
%   tjsr_zeroJsr(M)  decides if the JSR of a set of matrices M is equal to zero.
%
%   [OUT,NORMHIST] = tjsr_zeroJsr(M)
%                   For M a cell array, tjsr_zeroJsr starts with a basis of 
%                   R^N and determines if for all possible products of
%                   length N the image of X is equal to zero. The function
%                   assumes X=0 when norm(X) < 1e-6.
%
%   OUTPUTS:
%
%   - OUT       is a boolean, equal to true if the JSR of M is zero, and 
%               false if the JSR of M is > 0.
%
%   - NORMHIST  is an array of length N, which contains the norm of the
%               image of the basis X after each step.
%
%   [...] = JSR_ZEROJSR(M,TOL)
%                   The value TOL is used to decide if X=0. When
%                   norm(X)<TOL, the function assumes X=0.
%
%    R.Jungers,
%      "The Joint Spectral Radius: Theory and Applications"
%      Vol. 385 Springer Science & Business Media, 2009.
%
% Original version: Jungers, 2015
% Modified by: tommsch, February 2019

%#ok<*ALIGN>

if(nargin < 2)
    tol = 1e-6; end;

if(~isnumeric(tol) || tol < 0 || imag(tol)~=0)
    error('Input "tol" must be a real positive value.'); end;

m = length(M);
n = length(M{1});
out = false;

X = eye(n);
normHist = zeros(m);

for i=1:n
    temp = zeros(n);
    for j=1:m
        temp = temp + M{j}'*X*M{j}; end;
    X = temp;
    normHist(i) = norm(X);
    if(normHist(i)<tol) 
        out = true;
        return; end;
end

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 