function [C, idxminnew, idxmaxnew] = equalizeminidx(varargin);
% [C, idxminnew, idxmaxnew] = equalizeminidx({A1, idxmin1; A2, idxmin2; ...; An, idxminn}}, [options]);
% Makes a set of matrices the same size
%
% Input:
%   N x 2 cell array
%   first row: matrices, second row: index of first entry of matrix
%   Indices must be given as a column vector, all of their sizes must be equal
%
% Output:
%   Zero padded array such that index of first entry is for all matrices the same,
%   and such that their size is equal
%
% See also: setminidx, sequence
%
% E.g.: vdisp(equalizeminidx({[1 1 2; 1 2 3],[1 1 0]';[1],[0 0 1]'}))
%
% Written by: tommsch, 2018


A=varargin{1};


%find dimension
dim=size(A{1,2},1);
for i=1:size(A,1)
    dim=max(dim,size(A{i,2},1));
end

%find minimal minidx
MIN=ones(dim,1)*inf;
MAX=ones(dim,1)*-inf;
for i=1:size(A,1);    
    MIN=min(MIN,A{i,2});
    MAX=max(MAX,A{i,2}+sizem(A{i,1},[],dim).');
end;
MAX=MAX-1;

C=cell(size(A,1),1);
for i=1:size(A,1)
    C{i}=setidx(A{i,1}, A{i,2},MIN,MAX);
end

idxminnew=MIN;
idxmaxnew=MAX;

end