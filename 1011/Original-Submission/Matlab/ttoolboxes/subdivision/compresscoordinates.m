function [ xyzv ] = compresscoordinates( varargin )
% [ xyzv ] = compresscoordinates( c, M, [options] )
% Returns an array representing the graph of a function.
% i.e.: The points in c are at integer values
%       The points in xyzv are at M^{-1}(idx), or precisely
%       xyzv(:,i)= [M^{-1}*ind2subn(size(c),i); c(i)] %where i is the linear index of c
%
% Input:
%   c               (dim-array) Sequence with the function values
%   M               (dim x dim matrix) Matrix which defines where the function values are
%
% Options
%   'idx',val       Index of the first entry in c. Default minidx=0;
%   'verbose',val   Verbose level
%
% Output:
%   xyzv            (dim x N matrix) List of arguments/function values. Can be plotted by plotn
%
% E.g.: compresscoordinates([10 20 30; 40 50 60; 70 80 90],[2 1; -1 1])
%
% See also ind2subm, plotm
%
% Written by: tommsch, 2016

c=varargin{1};
M=varargin{2};
dim=size(M,2);
idx=parsem('idx',varargin,zeros(dim,1)); 
%removezeros=parsem('removezeros',varargin,zeros(dim,1)); 

IDX=ind2subm(sizem(c),1:numel(c))-1;
IDX=setplus(IDX,idx,'nounique');
xyzv=[M\IDX; c(:).'];

end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.