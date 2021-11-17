function S = multiplyS(varargin)
% [ S ] = multiplyS( S1, S2, ..., Sn, [options] )
% Concatenates subdivision operators, i.e. S = S1*S2*...*Sn.
% ! This function is experimental !
% 
% Input:
%   Si              Subdivision operators
%
% Options:
%   'verbose',val   Verbose level.
%   'dim',val       Dimension of subdivision operators. Only used if there are no arguments at all.
%   
% Output:
%   S       S1*S2
%
% Note:
%   If there are no arguments the function returns the idenditiy subdivision operator: a=1; M=eye(dim); D=0;
%
% E.g.: S1=getS('1_143');S2=getS('1_1133'); S=multiplyS(S1,S2);
% 
% Written by: tommsch, 2018

% XX There is something wrong with the ordering of multiple schemes.
% XX Output of [c,PM,x]=blf(multiplyS(S2,S1),'iteration',1);
% XX and       [c12,PM12,x12]=blf({[],[1 2]},[S2;S1],'iteration',2);
% XX is not the same

 %#ok<*ALIGN>

%Parse input
[verbose,varargin] = parsem({'verbose','v'},varargin,1);
[dim,varargin] = parsem('dim',varargin,1);
if(~all(cellfun(@isS,varargin))); 
    error('Input must be subdivision schemes'); end;
if(size(varargin,2)==1); 
    S=varargin{1}; 
    return; end;
if(size(varargin,2)==0); 
    S=getS({1,eye(dim),0}); 
    return; end;


S=multiplyS_worker(varargin{1}, varargin{2});
for i=3:size(varargin,2)
    S=multiplyS_worker(S,varargin{i});
    
end

if(rho(S{2}^(-1))>1); 
    warning( 'multiplyS:notexpanding', 'Output-scheme has non-expanding dilation matrix, i.e. is no subdivision scheme.' ); end;

end

function S = multiplyS_worker(S1, S2)

%a1=S1{1}; a2=S2{1};
M1=S1{2}; M2=S2{2};
D1=S1{3}; D2=S2{3};

a12=blf([1 2],[S1;S2],'plot',0,'verbose',0,'iteration',2);
M12=M2*M1;
D12=setplus(M2*D1,D2);

sze=size(S1,2);
S = getS({a12,M12,D12});
S{1,sze}='multiplied';

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 