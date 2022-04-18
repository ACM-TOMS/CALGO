function [ c, ic ] = smallestchoice( varargin )
% [ c, ic ] = findsmp( C, [options] )
% Tries to find an optimal choice of representatives, namely the following problem:
% Sought is minimal a set of elements such that from each set in C at least on element is chosen,
% and under all such choices c1,...,ck we sought for a choice such that max(v1,...,vk) is minimal.
%
% The algorithm used is greedy and does not guarantee the correct solution of the posed problem.
%
% Input:    
% ======
%   C                       1xN cell array of dx1 column vectors, 
%
% Options:
% ========
%   'equalfun',h            handle to function taking two vectors of length d, returning true/false, default: equalfun = @(x,y) isequal(x(1:end-1),y(1:end-1)); OR
%                           double, in which case equalfun = @(x,y) norm(x(1:end-1)-y(1:end-1))<double;
%                           Function which determines which ci are considered to be equal
%   'valfun',h              handle to function taking vector of length d, returning double, default: valfun = @(x) x(end);
%                           Function which determines the values vi.
%
% Output:
% =======
%   c                       1xN vector of chosen elements
%   ic                      cell array of linear indices  chosen elements
%
% Written by: tommsch, 2020

%row-indices
EL = 1; %elements %values of EL, VAL, IC are fixed and cannot be changed
VAL = 2; %value of element
IC = 3; %linear index
CUMVAL = 4; %valfunc of all same elements
SET = 5; %original set


%parse input
[equalfun,varargin] = parsem( {'equalfun','ef'}, varargin, [] ); %Default: @(x,y) norm(x(1:end-1)-y(1:end-1))<1e-8
[valfun,varargin]   = parsem( {'valfun','vf'}, varargin, [] ); %@(x) x(end)
[verbose,varargin]  = parsem( {'verbose','v'}, varargin, 1 ); %@(x) x(end)
Co = varargin{1}; 
varargin(1) = [];
parsem( varargin, 'test' );
if( isempty(Co) );
    c = {[]};
    ic = {[]};
    return; end;

%pre-process input
vprintf( 'Preprocess input.\n', 'imp',[2 verbose] );

n = sum(cellfun(@(x) size(x,2), Co));
C = nan( 1, n );
id = 1; 
for i = 1:n
    if( ~isnan(C(i)) );
        continue; end;
    C(i) = id;
    id = id+1;
    if( isa(equalfun,'function_handle') )
        val = num2cell( [Co{:}], 1 );
        idx = [zeros( 1, i ) cellfun( @(x) equalfun(val{:,i},x), val(i+1:end) )];
    elseif( isscalar(equalfun) )
        val = num2cell( [Co{:}], 1 );
        idx = [zeros( 1, i ) cellfun( @(x) norm(val{:,i}(1:end-1)-x(1:end-1))<equalfun, val(i+1:end) )];
    elseif( size(Co{1},1)>2 );
        val = [Co{:}];
        idx = [zeros( 1, i ) ismember( val(1:end-1,i+1:end).', val(1:end-1,i).', 'rows' ).']; 
    else
        C = [Co{:}];
        C = C(1,:);
        break;
        end; 
    C(logical( idx )) = C(i);
end;

if( ~isempty(valfun) );
    val = num2cell( [Co{:}], 1 );
    V = cellfun( valfun, val );    
else
    V = [Co{:}];
    V = V(end,:); end;

ic = cellfun( @(x) 1:size(x,2), Co, 'UniformOutput',false );
ic = [ic{:}];

C = [C;V;ic];
C = mat2cell( C, 3, cellfun(@(x)size(x,2),Co) );


vprintf( 'Choose elements input.\n', 'imp',[2 verbose] );
%choose of duplicates best element in each partition
for i = 1:numel( Co );
    C{i}(SET,:) = i; %original set
    val = [];
    while( ~isempty(C{i}) );
        idx_c = C{i}(EL,1)==C{i}(EL,:);
        idx_v = find( C{i}(VAL,:)==min(C{i}(VAL,idx_c)) & idx_c, 1 );
        val = [val C{i}(:,idx_v)]; %#ok<AGROW>
        C{i}(:,idx_c) = []; end; %this line is slow
    C{i} = val; end;

%make everything to one matrix
C = horzcat( C{:} );

%compute value of elements
C(CUMVAL,:) = nan;
for i = 1:size(C,2);
    if( ~isnan(C(CUMVAL,i)) );
        continue; end;
    idx = C(EL,i)==C(EL,:); %find equal elements
    C(CUMVAL,idx) = max( C(VAL,idx) ); end;

C_full = C;
%select elements
chosen = []; %the selected elements
while( ~isempty(C) )
    occ = sum( C(EL,:)==C(EL,:)', 1 ); %number of occurences %XX this line is slow if the matrix C1 is large, maybe use hist
    idxocc = occ==max( occ ); %find elements of maximal occurence of not yet selected sets
    idxmincumval = C(CUMVAL,:)==min( C(CUMVAL,idxocc) ); %find elements with minimal value
    idx = find( idxocc & idxmincumval, 1 );
    el = C(EL,idx); %element with minimal value
    chosen = [chosen el]; %#ok<AGROW>
    idxel = C(EL,:)==el; %find in which sets the element occurss
    idxset = any( C(SET,:)==C(SET,idxel)', 1 ); %find chosen sets
    C(:,idxset) = []; end; %remove all sets which are chosen

vprintf( 'Chosen elements after greedy algorithm:              %r\n', chosen, 'imp',[2 verbose] );

C = C_full;
%find elements which can be removed
for i=1:numel(chosen) 
    %check if chosen without C1(EL,i) is also a valid choice
    chosen_i = chosen([1:i-1 i+1:end]);
    chosen_i = unique( C( SET, anym(C(EL,:)==chosen_i',1)) );
    if( isequal(chosen_i,1:size(Co,2)) )
        chosen(i) = nan; %#ok<AGROW>
        end; end;
chosen(isnan( chosen )) = [];

vprintf( 'Chosen elements after removing superfluous elements: %r\n', chosen, 'imp',[2 verbose] );

%find indices of chosen values and translate back indices to original values
ic = cell( 1, numel(Co) );
c = cell( 1, numel(Co) );
for i = 1:numel(Co);
    ic{i} = C(IC,anym(C(EL,:)==chosen.',1) & C(SET,:)==i);
    c{i} = Co{i}(:,ic{i}); end;



end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.