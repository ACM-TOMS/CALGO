function [ varargout ] = sizem(A, dim, N)
% [ varargout ] = sizem( A, [dim], [N] )
% Consistent behaviour of size, with regards to multi-dimensional applications.
%
% Input:
%       dim                 If dim>0 is given, then sizem has the same behaviour as size
%       N                   defines the length of the output vector. If dim and N are given, dim is ignored.
%                           If the dimension of A is larger than N, the output vector still has length of the dimension of A
%       both not given      consistent behaviour of @size regarding multi-dimensional applications
%
% Output:
%       If there is more than one output argument, then sizem has the same behaviour as size
%
% Remarks: 
%       Empty stuff has sizem==[]
%
% E.g.: sizem([2 3])
%       sizem([2;3])
%       sizem([2;3],[],3)
%       sizem(zeros(2,0,1))
%       sizem([])
%           
% See also: size
%
% Written by tommsch, 2018

if(nargin==2 || nargout>1);
    if(nargin==2)
        [varargout{1:nargout}]=size(A,dim);
    else
        [varargout{1:nargout}]=size(A);
    end
else
    out=size(A);
    if(nargin==3)
        out(1,size(out,2)+1:N)=1;
        out=out(1:N);
    elseif(isequal(out,[0 0]));
        out=[];
    elseif(out(2)<=1 && size(out,2)==2); 
        out=out(1); 
    end;
    
    varargout{1}=out;
end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 