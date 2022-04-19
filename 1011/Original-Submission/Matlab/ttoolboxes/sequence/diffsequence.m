function [cout, mu] = diffsequence( varargin )
% [ cout, mu ] = diffsequence( cin , DIFF, options)
% Computes backward differences of sequences.
%
% Input:
%   cin       multidim. array. 
%   DIFF      number or vector
%
% Options:
%   'cell'          return value is always a cell array
%   'equalsize'     return arrays are zero-padded so that they have the same size
%   'dim'           the dimension of the input array. If not given, it is determined automatically
%
% Output:
%   c       If DIFF is a vector, then the DIFF^th backward difference is computed
%           If DIFF is a number, then all backward differences up to order DIFF are computed and returend as a cell array of equal size (zero padded)
%   mu      If DIFF is a vector, then the d is a cell array and contains the direction of the corresponding backward difference in c
%           If DIFF is a number, then d == DIFF 
% Note:
%   Be carefull if cin is a vector, the behaviour could be different than you expect. See the example below.
%   Since the first index of the sequences does not change, there is no need to deal with this index.
%
% E.g.: diffsequence([1 2 3; 4 5 6],1,'equalsize')
%       diffsequence([1 2 3]',1,'cell')
%       diffsequence([1 2 3]',1)
%       diffsequence([1 2 3],[0;1])
%
% See also: diff, sequence
%
% Written by: tommsch, 2018

% 2019-03: Added option 'dim'

cin=varargin{1};
DIFF=varargin{2};
dim=parsem('dim',varargin,ndimsm(cin));
cellflag=parsem('cell',varargin);
equalsize=parsem('equalsize',varargin);

DIFF=squeezem(DIFF);

if(~isscalar(DIFF)); %the direction for diffsequ
    if(size(DIFF,2)~=1); 
        error('diffsequ: DIFF must be a vector with dim entries or a scalar.'); end;
    mu=DIFF; %if vector, then DIFF is the direction
elseif(dim==1);
    mu=DIFF;
else;
    mu=constructmu(DIFF,dim);
end;

cout=cell(size(mu,2),1);
d=cell(size(mu,2),1);
for i=1:size(mu,2);
  cout{i}=diffsequ_worker(cin,mu(:,i)');      
  d{i}=mu(:,i)';
end;

if(equalsize);  %make cell arrays the same size
    MAX=max(cell2mat(cellfun(@size,cout,'UniformOutput',0)),[],1); %find biggest index of all cell arrays
    MAX=num2cell(MAX);
    for i=1:size(cout,1)
        if(anym(sizem(cout{i})<[MAX{:}]));
            cout{i}(MAX{:})=0;
        end
    end
end    

if(size(cout,1)==1 && ~cellflag);   %return matrix if cell array consists only of 1 element
    cout=cout{1}; 
end


end

function [c] = diffsequ_worker(c, DIFF)
    c=padarraym(c,DIFF);
    for i=1:size(DIFF,2); 
        if(DIFF(i)==0); 
            continue; end;
        c=diff(c,DIFF(i),i);     
    end 
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 