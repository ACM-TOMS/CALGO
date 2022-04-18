function [z] = ordering2num(varargin)
% [ z ] = ordering2num( oo, S, [options] )
% Takes an ordering and computes the corresponding number in \RR^dim.
%
% Input:
%   S                   subdivision operators
%   oo                  (ordering) sequence of digits to be chosen. 
%                       Format of oo: Either {2xN,2xM} : first row ordering of subdivision operators, second row ordering of digits,
%                                     or {1xN,1xM} : first row ordering of digits
%                       First cell is the non-periodic part, second cell is the periodic part
%
% Options: 
%   'verbose',val       Verbose level
%
% Output:
%   z                   column-vector, the corresponding number
%
% Info: This function is much faster than num2ordering. Thus if possible this function should be used.
%
% E.g.: ordering2num({[2 1 ],[]},getS({[],2}))
%       ordering2num({[2 2 ],[2 1 2]},getS({[],10}))
%
% See also: num2ordering, getS
%
% Written by: tommsch, 2018

% XX Wenn 0 in einem digit set nicht drinnen ist, und der periodische Teil leer ist, setze den periodiscen Teil auf die Äquivalenzklasse von 0.
% XX Gib trotzdem eine Warnung aus.

%#ok<*ALIGN>

verbose=    parsem({'verbose','v'},varargin,1);



%parse second input
if(isnumeric(varargin{2}));
    M{1}=varargin{2};
    if(isscalar(M{1}));  
        D{1}=0:M{1}-1; 
    else; 
        D{1}=constructdigit(M{1}); end;
    J=1;
    dim=size(M{1},2);
    ZERO=zeros(dim,1);
else    
    S = getS(varargin{2});
    J=size(S,1);
    dim=size(S{1,2},1);
    ZERO=zeros(dim,1);
    M=S(:,2)';
    D=S(:,3)';
end

%parse oo
oo = varargin{1};


if(J>1 && ( size(oo{1},1)==1 || size(oo{2},1)==1) );  %make checks
    error('ordering2num: ''do'' needs a 2xN array.'); end;

if(size(oo,2)~=2); 
    error('ordering2num: ''oo'' needs to be 2x1 cell array.'); end;
L={size(oo{1},2),size(oo{2},2)};
if(size(oo{1},1)==1); 
    oo{1}=[ones(1,L{1}); oo{1}]; end;
if(size(oo{2},1)==1); 
    oo{2}=[ones(1,L{2}); oo{2}]; end;



%compute periodic part first
z=ZERO;
if(~isempty(oo{2}))
    for i=size(oo{2},2):-1:1
        j=oo{2}(1,i); %which subdiv operator
        d=oo{2}(2,i); %which digit
        z=M{j}\(z+D{j}(:,d));
    end
    z=tbuildproduct(M,oo{2}(1,:))/(tbuildproduct(M,oo{2}(1,:))-eye(dim))*z; %periodic part
    if(~isempty(oo{1})); 
        z=tbuildproduct(M,oo{1}(1,:))\z; end; %shift periodic part
else
    for i=1:J; 
        if(isempty(intersect(D{i}',ZERO','rows')));
            vprintf('ordering2num: ''do{2}'' can''t be empty, since 0 is not a digit for all subdivision operators.\n','cpr','err','imp',[0 verbose]);
        end;
    end;
end

if(~isempty(oo{1}))    
    for i=size(oo{1},2):-1:1
        j=oo{1}(1,i); %which subdiv operator
        d=oo{1}(2,i); %which digit
        z=tbuildproduct(M,oo{1}(1,1:i))\D{j}(:,d)+z;
    end
end

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 


