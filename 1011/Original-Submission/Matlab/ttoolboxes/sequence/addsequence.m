function [C, cmin]=addsequence(A,B,amin,bmin)
% C = addsequence(A, B, [amin], [bmin])
% Sums matrices of unequals size
%
% Input: 
%   A,B         arrays, can be empty
%   amin,bmin   default=0, index of the first entry as column vector
%
% Output:
%   C           A+B
%   cmin        index of first entry of C as column vector
%
% E.g.: addsequence([1 2 3],[4;5;6])
%       addsequence([1 2 3],[4;5;6],[1;1],[1;-1])
%
% See also: sequence
%
% Written by: tommsch, 2017

dim=max(length(sizem(A)),length(sizem(B)));
if(nargin==2); 
    amin=zeros([dim 1]); bmin=amin;  %sic. the '1' ensures that I dont get a square matrix
end;

[C,cmin]=equalizeminidx({A,amin; B,bmin});
C=cat(3,C{:});
C=sum(C,3);

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 