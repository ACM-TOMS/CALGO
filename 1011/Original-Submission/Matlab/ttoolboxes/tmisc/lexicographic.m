function [ZOUT] = lexicographic(varargin)
% [ zout ] = lexicographic( zin, [norm] )
% Orders vectors in a norm-lexicographic ordering.
%
% Input:
%   zin             Array of column vectors. The columns are sorted first by their norm, second by their lexicographic-value
%   [norm]          either a number or 'inf'. The norm used to sort
%
% Output:
%   zout            The sorted array
%
% Example: lexicographic( [0 1 0; 0 1 2; 3 1 1] ,2)
%
% See also: reducelength
%
% Written by tommsch, 2016

 %#ok<*ALIGN>

ZZ=varargin{1};

if(nargin>=2);
    norm=varargin{2};
    if(isscalar(norm) && isnumeric(norm));                 
        NORMS=sum(abs(ZZ).^norm).^(1/norm);     %p-norm (p='norm')
    elseif(strcmp(norm,'inf'));         
        NORMS=max(abs(ZZ));                     %inf-norm
    else;                               
        error('Wrong value for ''norm''.'); end;
else;                               
    NORMS=sum(abs(ZZ),1);                     %1-norm (default)
end;

EPS=1e10;
NORMS=round( NORMS*1e10)/1e10;
l=1;
if(any(isnan(NORMS))); 
    error('Lexicographic Ordering failed. NaN''s occured.\n'); end;

for i=unique(NORMS);
    if(i==inf);         
        ZC=ZZ(:,NORMS==inf);
    else;               
        ZC=ZZ(:,abs(NORMS-i)<1/EPS);
    end;
    
    ZC=sortrows(ZC')';
    ZOUT(:,l:l+size(ZC,2)-1)=ZC;
    l=l+size(ZC,2);
end
if(size(ZOUT,2)~=size(ZZ,2)); 
    disp(ZZ);
    disp(ZOUT);    
    error('Lexicographic Ordering failed. Output Size \neq Input size.\n'); end;

end