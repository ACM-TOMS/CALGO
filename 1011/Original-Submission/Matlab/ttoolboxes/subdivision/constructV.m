function [V, Om] = constructV(varargin);
% [ V, Om ] = constructV( Om, [k] )
% Constructs a basis for the space V_k(Om), as described in Charina, Mejstrik, 2018.
%
% Input: 
%   Om          The set for which V_k shall be constructed as an array of column vectors
%   k           index k, optional. Either an integer or a vector of integers (greater equal zero)
%               If k is empty, then all spaces V_k, where dim V_k>0 are computed
%
% Options:
%   '01'        Give input as a 0/1 matrix
%   'verbose'   Verbose level
%
% Output:
%   V           matrix/or cell array of matrices, each column is one basis vector for the space V_i i \in k
%   Om          The set for which V_k is constructed as an array of column vectors.
%
% Note:

%
% E.g.: vdisp(constructVt([0 1 2 3]))
%       Om=constructOmega('2_butterfly'); constructV(Om,1)
%       constructV([1 0 1; 0 1 1; 1 0 1],1,'01')
%
% See also: constructVt, constructU
%
% Reference:
%   Set V_k is described in:
%       Maria Charina, Thomas Mejstrik,
%       Multiple multivariate subdivision schemes: Matrix and operator approaches, 
%       Journal of Computational and Applied Mathematics, 2018
%
% Written by: tommsch, 2016

%#ok<*ALIGN>

% Parse Input
Om=varargin{1};
if(isempty(Om));  
    if(nargin==1) 
        V={[]}; 
        return;
    else;
        V=[];
        return; end;
end;
    

[zeroone,varargin]=parsem('01',varargin);
[verbose,varargin]=parsem({'verbose','v'},varargin,1);
if(zeroone); 
    dim=length(size(Om));
    Om=supp(Om,dim,zeros(dim,1));
end;
vprintf('Om: \n%v\n',Om,'imp',[2 verbose]);

if(size(varargin,2)>=2 && issym(varargin{2}));
    k=varargin{2};
    varargin(1:2)=[];
    kflag=1; %flag if 'k' is given or not
elseif(size(varargin,2)>=2 && ~isempty(varargin{2}) && isnumeric(varargin{2}) && ~ischar(varargin{2})); %test if second argument is a valid number for k
    k=varargin{2};
    kflag=1;%flag if 'k' is given or not
    varargin(1:2)=[];
else;
    k=0:size(Om,2);
    kflag=0; %flag if 'k' is given or not
    if(size(varargin,2)>=2 && isempty(varargin{2})); varargin(2)=[]; end; 
    varargin(1)=[];
end
parsem(varargin,'test');

j=0; %idx
V=cell(1,size(k,2)); %allocate space big enough
while(true);
    j=j+1;
    V{j}=constructV_worker(Om,k(j));
    if( (kflag && j==size(k,2)) || (~kflag && isempty(V{j})) );        
        break; end; %test stopping criterion
end
if(~kflag); 
    V=V(~cellfun('isempty',V)); end;
if(size(k,2)==1); 
    V=V{1}; end;


end


function [ V ] = constructV_worker( Om, k)
%for Omega a matrix of vectors


if(k>8); vprintf('constructV: Function may return false values for k>8.\n','cpr',[.6 .4 0]); end;

sizeOm=size(Om,2);
dim=size(Om,1);
mu=mixvector(0:k,dim);  %all possible multiindices whose sum is less than k
[~,osum]=sort(sum(mu,1));
mu=mu(:,osum);
mu=mu(:,sum(mu,1)<=k & sum(mu,1)>=0);
mu=fliplr(mu);
MV=zeros(size(mu,2),sizeOm);
for j=1:size(mu,2);
    for i=1:sizeOm;
        MV(j,i)=prod(Om(:,i).^mu(:,j));
    end
end
V=null(MV,'r'); 


end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.

