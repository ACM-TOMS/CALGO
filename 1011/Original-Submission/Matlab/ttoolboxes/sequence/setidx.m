function [A, idxminnew, idxmaxnew] = setidx(varargin);
% Anew = setminidx(A, idxmin, idxminnew, [idxmaxnew]);
% Changes the upper left index of an array by padding the array.
%
% Input:
%   A                   the array to be padded
%   idxmin              (column vector) the index of the first entry in A
%   idxminnew           (column vector) the new index of the first entry in Anew. idxminnew<=idxmin
%   idxmaxnew           (optional, column vector) the new index of the last entry in Anew. idxmaxnew>="idxmaxold"
%
% Output:
%   A                   The padded array
%   idxminnew           Equals always idxmin
%   idxmaxnew           The idx of the last entry in A
% 
% E.g.: setidx([1 2; 2 3],[1 1]',[0 0]',[4 3]')
%
% See also: equalizeminidx, padarraym
%
% Written by: tommsch, 2017


A=varargin{1};

idxmin=varargin{2};
idxminnew=varargin{3};
if(~iscolumn(idxmin) || ~iscolumn(idxminnew)); 
    error('indices must be given as column vector.'); end;
dim=size(idxmin,1);

padmin=idxmin-idxminnew;
if(any(padmin<0)); 
    error('idxminnew>idxmin'); end;

A=padarraym(A,padmin.','pre'); %can handle cells and symbolic stuff

if(size(varargin,2)>=4);
    idxmax=sizem(A,[],dim).'+idxminnew-1;
    idxmaxnew=varargin{4};
    padmax=idxmaxnew-idxmax;
    if(any(padmax<0)); 
        error('idxmax>idxmaxnew'); end;
    A=padarraym(A,padmax.','post');
else
    idxmaxnew=idxminnew+sizem(A)-1;
end


end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 