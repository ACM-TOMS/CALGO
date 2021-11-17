function [ Vt, Om, Xmuf ] = constructVt(varargin);
% [ Vt, Om, Xmuf ] = constructVt( Om, [k], [options] )
% Constructs a basis for the space Vtilde_k(\Omega) as described in Charina, Mejstrik, 2018.
%
% Input: 
%   Om          The set for which Vtilde_k shall be constructed as an array of column vectors
%   k           index k, optional. Either an integer or a vector of integers (greater equal zero)
%               If k is empty, then all spaces Vtilde_k, where dim Vtilde_k>0 are computed
%
% Options:
%   '01'        give Omega as a 0/1 matrix
%
% Output:
%   Vt          matrix/cell array of matrices. Each column is one basis vector for the space Vtilde_i, i \in k
%   Om          The set for which Vt_k is constructed as an array of column vectors.
%   Xmuf        scalar/vector. If set, all sets X_mu are non-empty, see [Mejstrik, PhD-thesis, 2019].
%
% E.g.: vdisp(constructVt([0 1 2 3]))
%       Om=constructOmega('2_butterfly'); constructVt(Om,1)
%       constructVt([1 0 1; 0 1 1; 1 0 1],1,'01')
%
% See also: constructV, constructU
%
% Reference:
%   Set Vtilde_k is described in:
%       Maria Charina, Thomas Mejstrik,
%       Multiple multivariate subdivision schemes: Matrix and operator approaches, 
%       Journal of Computational and Applied Mathematics, 2018
%
% Written by: tommsch, 2016

% XX written very ugly, could be rewritten
% XX Om=[     0     0     0     1     1     1     2     2     2;     0     1     2     0     1     2     0     1     2];
%    constructVt(4); %makes 0x9 array instead of 9x0 array

%#ok<*ALIGN>

Om=varargin{1};
if(isempty(Om));  
    if(nargin==1) 
        Vt={[]}; 
        return;
    else;
        Vt=[];
        return; end;
end;
[zeroone,varargin]=parsem('01',varargin);
[verbose,varargin]=parsem({'verbose','v'},varargin);

if(zeroone); 
    dim=length(size(Om));
    Om=supp(Om,dim,zeros(dim,1));
end;
vprintf('Om=\n%v\n',Om,'imp',[2 verbose]);

if(size(varargin,2)>=2 && ~isempty(varargin{2}) && isnumeric(varargin{2}) && ~ischar(varargin{2})); %test if second argument is a valid number for k
    k=varargin{2};
    kflag=1;%flag if 'k' is given or not
    varargin(1:2)=[];
else;
    k=0:size(Om,2);
    kflag=0;%flag if 'k' is given or not
    varargin(1)=[];
end
parsem(varargin,'test');

j=0; %idx
Vt=cell(1,size(k,2)); %allocate space big enough
Xmuf=zeros(size(Vt)); %allocate space big enough

while(true);
    j=j+1;
    [Vt{j},Xmuf(j)]=constructVt_worker(Om,k(j));
    if( (kflag && j==size(k,2)) || (~kflag && isempty(Vt{j})) );        
        break; end; %test stopping criterion
end
if(~kflag); 
    Vt=Vt(~cellfun(@isempty,Vt)); 
    Xmuf=Xmuf(1:numel(Vt)); end;
if(size(k,2)==1); 
    Vt=Vt{1}; end;


end

function [ Vt, Xmuf ] = constructVt_worker( Om_idx, k )
dim=size(Om_idx,1);
[Om,~,idx]=characteristic(Om_idx); %get characteristic function and sort Om
Vt=zeros(0,numel(Om)); % is getting transposed at the end

%construct all difference sequences
cc=diffsequence([1],k+1,'dim',dim,'cell');
Xmuf=zeros(1,numel(cc));

for i=1:numel(cc); %iterate over all mu's
    coriginal=cc{i};
    sizecoriginal=size(coriginal);
    sizeOmega=size(Om);
    if(length(sizecoriginal)<length(sizeOmega));
        sizecoriginal(length(sizeOmega))=1;
    end;
    sizecoriginal(sizecoriginal==0)=1;
    Loriginal=sizeOmega-sizecoriginal;
    if(any(Loriginal<0)); 
        continue; end;
    L=num2cell(Loriginal);
    for j=1:size(L,2);
        L{j}=0:L{j};
    end;
    Lstart=mixvector(L);
    Lend=setplus(Loriginal',-Lstart,'stable');

    for j=1:size(Lstart,2);
        c=coriginal;
        c=padarraym(c,Lstart(:,j)','pre');
        c=padarraym(c,Lend(:,j)','post');

        if (isequal((Om~=0).*(c~=0),c~=0));
            Xmuf(i)=1;
            Vt=[Vt; reshape(c,numel(c),1)']; %#ok<AGROW>
        end;
    end

end
Xmuf=all(Xmuf);
if(isempty(Vt)); 
    return; end;

Vt(:,find(~Om))=[]; %#ok<FNDSB>

RANK=rank(Vt);
for i=1:size(Vt,1);
    if(rank(Vt([1:i-1 i+1:end],:))==RANK); 
        Vt(i,:)=0;  
    end;
end   
Vt=removezero(Vt,1);
Vt=Vt';
Vt=Vt(idx,:);
end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.