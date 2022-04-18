function [ A, maxindex ] = mixvector( Z,dim,n )
% [ A, maxindex ] = mixvector( Z, dim, [n] )
% Constructs all possible combinations of values.
%
% Input: 
%   Z           values of which combinations are constructed, 
%                   vector
%                   (experimental) cell array
%   dim         length of combinations
%   n           optional. If not given returns all possible combinations. If given returns only the n^th combination
%
% Output:
%   A           the combinations
%   maxindes    number of possible combinations wrt. elements in Z and dim
%
% E.g.: mixvector(0:1, 4)
%       mixvector({[2 0],[3 0],[4 4]},3)
%
% Written by tommsch, 2017
% Uses the function allcomb from (c) Jos van der Geest, email: samelinoa@gmail.com

%#ok<*ALIGN>

if(iscell(Z)==1); 
    A=mixvector_cell(Z); 
    return; end
if(dim==1); 
    maxindex=size(Z,2);
    if(nargin==3 && n==0); 
        A=size(Z,2); 
        return; end;
    if(nargin==3 && n>maxindex); 
          error('Index exceeds array bounds.'); end;
    if(nargin==3); 
        A=Z(n); 
        return; end;
    A=Z;    
    return; end;
maxindex=length(Z)^dim;
if(nargin==3 && n==0); 
    A=maxindex; 
    return; end;

  if(nargin>=3);
      if(n>maxindex); 
          error('value of n too big'); end;
      number=dec2basem(n-1,length(Z),dim);
      A=Z(number+1)';
      return; end;

try; 
    MEM=tavailable_memory/2;
catch; 
    MEM=2*1024*1024*1024; %4GB of Data
end
if(dim*maxindex>MEM) 
    cprintf('blue','Do not use mixvector with big arguments.\n Your machine will slow down.\n');
    fprintf('Use instead 3 parameters: mixvector( Z,dim,n)=''mixvector(Z,dim)(:,n)''.\n');
    fprintf('You might want to cancel the operation now.\nPress any button to continue.');
    pause;
end;
A=zeros(dim,maxindex);
%s=whos('A');
%if(s.bytes>4*1024) %4KB of Data



for i=0:dim-1
    ROW=repmat(Z,length(Z)^i,1);
    ROW=reshape(ROW,1,[]);
    ROW=repmat(ROW,1,length(Z)^(dim-i-1));
    A(dim-i,:)=ROW;
end
if(nargin>=3);
    A=A(:,n); end; %just a fallback, if i comment out the upper stuff

end

function [ A ] = mixvector_cell( Z)
dim=size(Z,2);
if dim==1; A=Z{1}; 
    return; end

A=allcomb(Z{:})';
end

function A = allcomb(varargin)

% ALLCOMB - All combinations
%    B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
%    in the arrays A1, A2, ..., and AN. B is P-by-N matrix is which P is the product
%    of the number of elements of the N inputs. This functionality is also
%    known as the Cartesian Product. The arguments can be numerical and/or
%    characters, or they can be cell arrays.
%
%    Examples:
%       allcomb([1 3 5],[-3 8],[0 1]) % numerical input:
%       % -> [ 1  -3   0
%       %      1  -3   1
%       %      1   8   0
%       %        ...
%       %      5  -3   1
%       %      5   8   1 ] ; % a 12-by-3 array
%
%       allcomb('abc','XY') % character arrays
%       % -> [ aX ; aY ; bX ; bY ; cX ; cY] % a 6-by-2 character array
%
%       allcomb('xy',[65 66]) % a combination
%       % -> ['xA' ; 'xB' ; 'yA' ; 'yB'] % a 4-by-2 character array
%
%       allcomb({'hello','Bye'},{'Joe', 10:12},{99999 []}) % all cell arrays
%       % -> {  'hello'  'Joe'        [99999]
%       %       'hello'  'Joe'             []
%       %       'hello'  [1x3 double] [99999]
%       %       'hello'  [1x3 double]      []
%       %       'Bye'    'Joe'        [99999]
%       %       'Bye'    'Joe'             []
%       %       'Bye'    [1x3 double] [99999]
%       %       'Bye'    [1x3 double]      [] } ; % a 8-by-3 cell array
%
%    ALLCOMB(..., 'matlab') causes the first column to change fastest which
%    is consistent with matlab indexing. Example: 
%      allcomb(1:2,3:4,5:6,'matlab') 
%      % -> [ 1 3 5 ; 1 4 5 ; 1 3 6 ; ... ; 2 4 6 ]
%
%    If one of the arguments is empty, ALLCOMB returns a 0-by-N empty array.
%    
%    See also NCHOOSEK, PERMS, NDGRID
%         and NCHOOSE, COMBN, KTHCOMBN (Matlab Central FEX)

% Tested in Matlab R2015a
% version 4.1 (feb 2016)
% (c) Jos van der Geest
% email: samelinoa@gmail.com

% History
% 1.1 (feb 2006), removed minor bug when entering empty cell arrays;
%     added option to let the first input run fastest (suggestion by JD)
% 1.2 (jan 2010), using ii as an index on the left-hand for the multiple
%     output by NDGRID. Thanks to Jan Simon, for showing this little trick
% 2.0 (dec 2010). Bruno Luong convinced me that an empty input should
% return an empty output.
% 2.1 (feb 2011). A cell as input argument caused the check on the last
%      argument (specifying the order) to crash.
% 2.2 (jan 2012). removed a superfluous line of code (ischar(..))
% 3.0 (may 2012) removed check for doubles so character arrays are accepted
% 4.0 (feb 2014) added support for cell arrays
% 4.1 (feb 2016) fixed error for cell array input with last argument being
%     'matlab'. Thanks to Richard for pointing this out.

narginchk(1,Inf) ;

NC = nargin ;

% check if we should flip the order
if(ischar(varargin{end}) && (strcmpi(varargin{end},'matlab') || strcmpi(varargin{end},'john')));
    % based on a suggestion by JD on the FEX
    NC = NC-1 ;
    ii = 1:NC ; % now first argument will change fastest
else
    % default: enter arguments backwards, so last one (AN) is changing fastest
    ii = NC:-1:1 ;
end

args = varargin(1:NC) ;
% check for empty inputs
if(any(cellfun('isempty',args)));
    warning('ALLCOMB:EmptyInput','One of more empty inputs result in an empty output.') ;
    A = zeros(0,NC) ;
elseif(NC > 1)
    isCellInput = cellfun(@iscell,args) ;
    if(any(isCellInput))
        if(~all(isCellInput))
            error('ALLCOMB:InvalidCellInput', 'For cell input, all arguments should be cell arrays.'); end;
        % for cell input, we use to indices to get all combinations
        ix = cellfun(@(c) 1:numel(c), args,'un',0) ;
        
        % flip using ii if last column is changing fastest
        [ix{ii}] = ndgrid(ix{ii}) ;
        
        A = cell(numel(ix{1}),NC) ; % pre-allocate the output
        for k=1:NC;
            % combine
            A(:,k) = reshape(args{k}(ix{k}),[],1) ;
        end
    else
        % non-cell input, assuming all numerical values or strings
        % flip using ii if last column is changing fastest
        [A{ii}] = ndgrid(args{ii}) ;
        % concatenate
        A = reshape(cat(NC+1,A{:}),[],NC) ;
    end
elseif(NC==1);
    A = args{1}(:) ; % nothing to combine

else % NC==0, there was only the 'matlab' flag argument
    A = zeros(0,0) ; % nothing
end

end


