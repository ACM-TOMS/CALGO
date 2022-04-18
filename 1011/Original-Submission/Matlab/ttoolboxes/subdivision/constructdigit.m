function [D]=constructdigit(varargin)
% [ D ] = constructdigit(M, [options])
% Constructs the usual digit set M[0,1)^s\cap \ZZ^s. 
%
% Input:
%       M                       a dilation matrix. I.e. a square matrix with all its eigenvalues greater than 1 in modulus.
%
% Options: 
%       'ZZ',val                If given, then only values of of the set val are used to compute the digit-set.
%                               The set val must be big enough, otherwise the output may not be a digit set.
%                               The function may restart itself, if val is too small.
%       'sym'                   the digits are computed symbolically, which takes a lot of time.
%                               If not given, but the output is not a digit set, then the algorithm may restart itself using the 'sym' option.
%       'random',val            default=0, returns a random digit set, with digits randomized at maximum by val
%       'classify'              returns a (size(ZZ,2) x abs(det(M)))-matrix. 
%                               'D' is the digit class where the entries from ZZ belong to. 
%       'verbose',val           verbose level
%
% Output:
%       D                       the digit set, each column is one digit
%
% Note:
%       Function has undefined behaviour if M is not an integer matrix
%
% E.g.: constructdigit( [2 0; 0 2] )
%
% Written by: tommsch, 2018

 %#ok<*ALIGN>

[ZZ,varargin]=         parsem({'ZZ','Z'},varargin,[]);
[verbose,varargin]=    parsem({'verbose','v'},varargin,1);
[symflag,varargin]=    parsem('sym',varargin);
[classify,varargin]=   parsem('classify',varargin); 
if(isempty(ZZ) && classify); 
    error('Option ''classify'' works only together with option ''ZZ'',<val>.'); end;
[randomval,varargin]=  parsem({'random','rand','r'},varargin,0); 
if(randomval && classify); 
    error('Option ''classify'' does not work together with option ''random''.'); end;
M=varargin{1}; 
varargin(1)=[];
parsem(varargin,'test');

dim=size(M,1);
m=int32(round(abs(det(M))));

if(m<1); 
    D=[]; 
    return; end;

if(m<2);
    RVAL=9;
    if(randomval); 
        D=randi(2*RVAL,dim,1)-RVAL; return;
    else; 
        D=zeros(dim,1);  return; end;
end;
    
 
if(symflag); 
    IM=inv(sym(M));
    epsilon=0;
else
    IM=inv(M);
    epsilon=abs(IM-round(IM))/2;
    epsilon=min(epsilon(epsilon>2*eps));
    
end;

D=zeros(dim,0);
if(~isempty(ZZ));
   maxindex=size(ZZ,2); 
   
   for i=1:maxindex
        v=ZZ(:,i);   
        dnew=v-M*floor(IM*v+epsilon/dim);
        D=[D dnew]; %#ok<AGROW>
   end   
else
    n=max(max(abs(M)));
    %n=norm(M,1);
    maxindex=mixvector(-n:n,dim,0);
    
    for i=1:maxindex
        v=mixvector(-n:n,dim,i); 
        dnew=v-M*floor(IM*v+epsilon/dim);
        D=[D dnew]; %#ok<AGROW>
    end
end

    

     
if(~classify)
    D=int32(D);
    D=unique(D','rows')';
    if(size(D,2)~=m && ~symflag); 
        vprintf('constructdigit: Fallback to symbolic computation - may need a long time.\n','cpr','err','imp',[ 2 verbose]);
        [D]=constructdigit(M, varargin{:},'sym'); return; end;

    D=double(D);
    if(randomval); 
        D=M*randi(randomval,dim,size(D,2))+D;
        D=D-round(summ(D)/numel(D));
    end;
else
    D=int32(D);
    if(~symflag); 
        vprintf('Classify may produce wrong output, if ''sym'' is not set.\n','cpr','err','imp',[ 1 verbose]); end;
end

end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.
