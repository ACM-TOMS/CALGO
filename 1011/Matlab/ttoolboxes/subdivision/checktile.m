function [ flag, C, T ] = checktile( varargin)
% [ flag, C, T ] = checktile( M ,D, [options])
% [ flag, C, T ] = checktile( S, [options] )
% Computes if a pair dilation-matrix/digit-set generates a tile (i.e. an attractor with area 1).
%
% Input:
%   M               dilation matrix
%   D               digit set for M
%
% Options:
%   'verbose',val   Verbose level
%   'legacy',val    Use old algorithm to compute the contact-vectors T, with set of integers of size [-val val]^dim.
%                   If val==0, a value for val will be determined automatically.
%   'basis',val     A basis for the lattice ZZ^2. If not given, the standardbasis is taken.
%
% Output:
%   flag            true:1, false:0, indefinite:NaN
%                       If flag==nan, the spectral radii of the contact matrix C should be computed by hand.
%   C               Contact matrix
%   T               Contact vectors
%
% Info: Lattice = \ZZ^dim
%       If D is not a digit set, the behaviour is undefined.
%
% E.g.: M1=[1 2; -2 -2];  D1=[0 1; 0 -1]; 
%       M2=[-1 1; -1 -2]; D2=[0 0 0; -2 -1 0]; 
%       M12=M2*M1; D12=setplus(M2*D1,D2);
%       C1=checktile(M1,D1);
%       C2=checktile(M2,D1);
%       C12=checktile(M12,D12);
%
% References:
% K.H. Groechenig and A. Haas, Self-Similar Lattice Tilings
%
% See also: tilearea
%
% Written by: tommsch, 2018

% Changelog: tommsch, 2019-05-08    Added option 'basis' to pass a basis for ZZ^2
%                                   Changed behaviour for verbose-level >= 2

% XX checktile(getS({[],[0 1;2 0],[0 0;1 0]'})); %yields wrong output. This is not a tile, but algorithm says so.
%tile(getS({[],[0 1;2 0],[0 0;1 0]'}))

%Parse input
%vprintf("Warning: Either the paper Groechenig/Haas contains an error, or there is a bug in the program.\nchecktile([0 1;2 0],[0 0;1 0]'); returns the wrong output.\n\n",'cpr','err');

[verbose,varargin]=parsem({'verbose','v'},varargin,1);
[legacy,varargin]=parsem('legacy',varargin,[]);
[basis,varargin]=parsem('basis',varargin,[]);

if(isS(varargin{1})); 
    M=varargin{1}{2};
    D=varargin{1}{3};
    varargin(1)=[];
else
    M=varargin{1};
    D=varargin{2};
    varargin(1:2)=[];
end

parsem(varargin,'test');

%Make contact vectors
if(isempty(legacy)); 
    T=contactvectors( M, D, basis ); 
else
    vprintf('Use old algorithm.\n,','imp',[2 verbose])
    T=contactvectors_old( M ,D, legacy, verbose ); 
end

vprintf('Digits: \n%v\n',D,'imp',[2 verbose])

vprintf('Contact vectors: \n%v\n',T,'imp',[2 verbose])

C=contactmatrix(M, D, T);
vprintf('Contact matrix: \n%v\n',C,'imp',[2 verbose])
if(verbose>=2);
    subplot(1,2,1);
    image(C,'CDataMapping','scaled'); 
    title('Contact Matrix'); 
    set(gca,'XAxisLocation','bottom','YAxisLocation','left');
    cbh=colorbar;
    set(cbh,'YTick',cbh.Limits(1):1:cbh.Limits(2))
    
    subplot(1,2,2)
    tile(getS({[],M,D}));
    title('Attractor');
end;


rhoval=max(abs(eig(C)));
detval=abs(det(M));
vprintf('Spectral radius of contact matrix= %i\nModulus of determinant of dilation matrix = %i\n',rhoval,detval,'imp',[2 verbose])

if(rhoval<detval-1000*eps*size(M,2)); %sure a tile
    vprintf('(M,D) generates a tile.\n','imp',[1 verbose],'cpr',[0.1 .5 .0]);
    flag=1;
elseif(norm(rhoval-detval)<.1)
    vprintf('Check spectral radius symbolically. This may take a while.\n','imp',[1 verbose]);
    rhoval=max(abs(eig(sym(C))));
    detval=abs(det(sym(M)));
    if(rhoval<detval)
         vprintf('(M,D) generates a tile.\n','imp',[1 verbose],'cpr',[0.1 .5 .0]);
        flag=1;
    elseif(rhoval>=detval)
        vprintf('(M,D) generates NO tile.\n','imp',[1 verbose],'cpr','err');
        flag=0;
    else
        vprintf('rho = %i \t(Spectral radius of contact matrix)\ndet = %i \t(Modulus of determinant of dilation matrix)\n',rhoval,detval,'imp',[1 verbose],'cpr','err');
        vprintf('(M,D) tiles if rho<det. The result is numerically indefinite.\n','imp',[1 verbose],'cpr','err');
        flag=NaN;        
    end
else
    vprintf('(M,D) generates NO tile.\n','imp',[1 verbose],'cpr','err');
    flag=0;
end



end

function [ C ] = contactmatrix( M, D, T );
nT=size(T,2);
C=zeros(nT,nT);
for k=1:nT
    for l=1:nT
        val=intersect( (M*T(:,k)+D).', (T(:,l)+D).', 'rows').';
        val=size(val,2);
        if(val); 
            C(k,l)=val;
        end;
    end
end
end

function [ T ] = contactvectors( M, D, basis )
%Computes the matrices $\mathcal{T}_n$ and $\mathcal{T}$ as described in 
%Algortihm from KH Groechnig and A Haas, Self-Similar Lattice Tilings, 
%Thm 2.2 and Lemma 4.5
dim=size(M,1);

%Make T0 set
if(isempty(basis))
    T0=unique([eye(dim) -eye(dim)].','rows').'; %start-vector for lattice = \ZZ^2, T0=\mathcal{T}_0  %I use unique to sort the vector
else
    T0=[basis -basis];
end


Tn{1}=T0;
detM=abs(det(M));

while(true)
    %vprintf('%v\n',Tn{end});
    val=M\setplus(Tn{end},D,-D);
    idx=allm(abs(val-round(val))<1/(2*detM),1);
    %Tn{end+1}=unique([round(val(:,idx)) Tn{end}].','rows').'; %#ok<AGROW> %In the paper it is only written, the Tn stabilize. They are not necessarily contained in each other!
    Tn{end+1}=unique([round(val(:,idx))].','rows').'; %#ok<AGROW> 
    if isequal(unique([Tn{1:end-1}].','rows'),unique([Tn{:}].','rows')); %#ok<ALIGN> %In the paper it is only written, the Tn stabilize. They are not necessarily contained in each other!
        T=unique([Tn{:}].','rows').';
        break; end;
end
idx=ismember(T.',zerosm(1,dim),'rows');
T(:,idx)=[]; %remove zero vector %See Definition of \mathcal{T}_n^* before Lemma 4.5
    

end


function [ T ] = contactvectors_old( M ,D, Zval, verbose ) 
%Computes the matrices $\mathcal{T}_n$ and $\mathcal{T}$ as described in 
%Algortihm from KH Groechnig and A Haas, Self-Similar Lattice Tilings, 
%Thm 2.2 and Lemma 4.5
dim=size(M,1);
%INDEX=[];

%Make Z lattice without origin
if(Zval==0); 
    Zval=(maxm(D(:))+1)*(max(abs(M(:)))+1)*size(D,1); 
    vprintf('Zval=%d\n',Zval,'imp',[1 verbose]);
end;
    
ZZ=mixvector(-Zval: Zval,dim); %represents the integers
zeroindex = ismember(ZZ.',zeros(1,dim),'rows');
ZZ=ZZ(:,~zeroindex);%remove (0,0)

%Make T0 set
T0=unique([eye(dim) -eye(dim)].','rows').'; %start-vector for lattice = \ZZ^2, T0=\mathcal{T}_0  %I use unique to sort the vector


Tn{1}=T0;

while(true)
    INDEX=[]; 
    for i=1:size(ZZ,2);
        INTERSECT=intersect(setplus(M*ZZ(:,i),D).',setplus(Tn{end},D).','rows');
        if(~isempty(INTERSECT));
            %ZZ(:,i)
            INDEX=[INDEX i]; %#ok<AGROW>
        end
    end
    Tn{end+1}=unique(ZZ(:,INDEX).','rows').'; %#ok<AGROW> %I use unique to sort the vector
    %Tn{end+1}=unique([Tn{end}.'; ZZ(:,INDEX).'],'rows').'; %#ok<AGROW> 
    if isequal(unique([Tn{1:end-1}].','rows'),unique([Tn{:}].','rows')); %In the paper it is only written, the Tn stabilize. They are not necessarily contained in each other! 
        T=unique([Tn{:}].','rows').'; break;
    end;    
end
    

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 

