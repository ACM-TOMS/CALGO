function I=intersectinterval(varargin);
% K = intersectinterval2(I_1, I_2, ..., [options]);
% Intersects intervals.
%
% Input: 
%   I_i,I2   (row-vector of two numbers OR a scalar OR empty array) the intervals, first entry must be smaller than second entry
%
% Options:
%   'failsafe',val      If intersection would be empty, val determines how to handle this case.
%                           0               no failsafe
%                           'minmax'        I(1)=min, I(2)=max
%                           'maxmin'        I(1)=max(min), I(2)=min(max)
%
% Output:
%   I       (row-vector of two numbers OR a scalar OR empty array) the intersection, I1 \cap I2
%
% Note: Input must not contain NaN.
%
% Eg: intersectinterval([4 5],[1 4],[-inf 4])
%
% See also: blockjsr
%
% Written by: tommsch, 2018

[failsafe,varargin]=parsem('failsafe',varargin,0);

if(size(varargin,2)==0); I=[-inf inf]; 
    return; end;
MIN=max(cellfun(@min,varargin));
MAX=min(cellfun(@max,varargin));

if(MIN>MAX);
    I=[];
elseif(searchincellarray(MAX,varargin,1) || searchincellarray(MIN,varargin,1)); 
    I=MAX; 
else;
    I=[MIN MAX];
end;

if(isempty(I));
    switch failsafe
        case 0; 
            %do nothing
        case 'maxmin'; 
            I(1)=MAX; I(2)=MIN;
        case 'minmax'; 
            I(1)=min(cellfun(@min,varargin)); I(2)=max(cellfun(@max,varargin)); 
        otherwise; 
            error('Unkown failsafe option');
    end;
end;
            
    

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   