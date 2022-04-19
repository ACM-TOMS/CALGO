function [ Cnest ] = unflatten(varargin)
% Cnest = unflatten(Cflat, C, [options])
% Changes a flat cell array to a nested cell array with topology given by C
% XX This function is slow, since it treats every cell entry by its own. This function should be rewritten.
%
% If this function throws an error, it could be, that the option 'assign' is missing.
%
% Input: 
%       C           a nested cell array which defines the topology of Cnest
%                   If C is not a cell array, then Cnest is not a cell array.
%       Cflat       a flat cell vector which contains the values of Cflat. Cflat and C must have the same number of elements in total.
%                   Otherwise, the behaviour is undefined.
%
% Options:
%   'assign'        All values of Cnest equal to Cflat 
%
% Output:
%       Cnest       the nested flat cell array, defined by the topoloy of Cflat and the values of C
%
% E.g.: C={[2 3 4] [5 ; 5 ; 5 ] {3 7; 4 6}; 40 50 {60  80}}
%       Cf=flatten(C)
%       D=unflatten({1 2 3 4 5 6 7 8 9 10},C)
%       vdisp(D)
%    
%       unflatten({10},[1])
%       unflatten(10,{1 2 3},'assign')
% 
% See also: flatten, nestedcellfun
%
% Written by: tommsch, 2017

%#ok<*ALIGN>

Cflat=varargin{1};
C=varargin{2};
assignflag=parsem('assign', varargin);

if(~iscell(C)) 
    if(iscell(Cflat)); 
        Cnest=Cflat{1}; 
    else
        Cnest=Cflat; end;
    return; end;

Cnest=C;
i=1; %counter for cells in Cdeep
j=1; %counter for elements in Cflat
while(true)
    if(isempty(C))
        %do nothing
    elseif(iscell(C{i}));
        if(~assignflag); 
            Cnest{i}=unflatten(Cflat(j:end),C{i});
        else
            Cnest{i}=unflatten(Cflat,C{i}); end;
        j=j+numel(C{i}); %update counter 
    else
        if(~assignflag); 
            Cnest{i}=Cflat{i};
        else
            Cnest{i}=Cflat; end;
        j=j+1;
    end
    i=i+1;
    if(i>numel(C)); 
        break; end;
end


end