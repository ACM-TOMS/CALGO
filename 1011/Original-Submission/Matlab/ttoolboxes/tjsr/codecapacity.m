function S=codecapacity(varargin)
% S=codecapacity(C, [options])
% Given forbidden difference patterns, returns a finite set of matrices, whose JSR is connected to its capacity.
% Input: 
%   C                   cell array of forbidden codes, each cell is a row vector of numbers from -base+1 to base-1
%                       all codes must have the same length and must be mutually different
%
% Options:
%   'base',val          (experimental) If not given, then val is computed by the numbers in the given code
%   'verbose',val      	verbose level
%   'plot',val          (default: false) Enables plot output            
%   
% Output:
%   Cell array of matrices, whose joint spectral radius shall be computed.
%
% E.g.: codecapacity({[1 1 0],[1 -1 0]})
%
% Written by: tommsch, 2018

verbose=parsem({'verbose','v'},varargin,1);
plotflag=parsem('plot',varargin,0);

C=varargin{1};
l=cellfun(@length,C);
base=parsem('base', varargin, max(cellfun(@(x) max(abs(x)),C)+1));
vprintf('Base: %i. ',base,'imp',[1 verbose]);
if(any(diff(l))); error('All codes must have the same length.'); end; %test if all C_i the same length
if(base~=2); vprintf('\nAlgorithm is only tested for base 2. Use algorithm with care.\n','cpr','err','imp',[1 verbose]); end;
if(base==2); vprintf('Maximal size of vertex cover to be computed: %i. ',base^sum(cellfun(@(x) base^nnz(x==0),C)),'imp',[2 verbose]); end;
%test if all codes are mutually different
l=l(1);
if(l<2); vprintf('\nLength of differences must be at least 2.\n','cpr','err','imp',[1 verbose]); end;
T=mixvector(0:base-1,l);
sze=size(T);

vprintf('Generate table of differences. ','imp',[1 verbose]);
D=reshape(T',[1 sze(2) sze(1)])-reshape(T',[sze(2) 1 sze(1)]); %table of differences
i1=cell(size(C)); i2=cell(size(C));
for i=1:length(C);
    C{i}=reshape(C{i},[1 1 length(C{i})]);
    [i1{i},i2{i}]=find(all(bsxfun(@eq, D, C{i}), 3)); %search for C{i} in the table of differences
end
D=[]; T=[]; %#ok<NASGU>
G=digraph(vertcat(i1{:}),vertcat(i2{:}));
i1=[]; i2=[]; %#ok<NASGU>

%make bipartite graph adjancy matrix
vprintf('Make bipartite graph adjancy matrix. ','imp',[1 verbose]);
P=zeros(base,base^(l-2));
P(1:base,1)=1;
A=P;
for i=1:base^(l-2)-1;
    A=[A; circshift(P,i,2)]; %#ok<AGROW>
end
A=repmat(A,[1 base]);
if(plotflag); subplot(2,1,1); plot(digraph(A),'Layout','subspace3');  title('Bipartite graph adjancy (de Brujin matrix)'); drawnow; end;

vprintf('Compute all locally minimal vertex covers. ','imp',[1 verbose]);
if(plotflag); subplot(2,1,2); plot(G,'Layout','force'); title('Graph for which vertex cover must be computed'); drawnow; end;
M=grVerCover(G);

vprintf('Generate matrices. \n','imp',[1 verbose]);
S=cell(size(M));
idx=find(A);
for i=1:numel(S)
    S{i}=A;
    S{i}(idx(M{i}))=0;
end

end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 
