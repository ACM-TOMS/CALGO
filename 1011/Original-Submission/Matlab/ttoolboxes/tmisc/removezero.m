function [ cout, origin ] = removezero( varargin)
% [c, origin] = removezero(c, what, ['keepdim')
% Deletes zeros in arrays in various ways.
%
% Input:  c                                         The array where zeros shall be removed
%
%         what  
%           'compact',row-vector of integers        Removes all zeros in the given directions, shifts non-zero elements to the beginning, and fills up with zeros afterwards, should work for nd-arrays
%           row-vector of integers                  Removes empty things iterativly in the directions given by the vector 
%           'border'                                Removes all zeros at the boundary, works for nd-arrays        
%           'outside'                               Removes all zeros at the boundary in planes which do not go through zero, works for nd-arrays
%           'inside'                                Removes all zeros at the boundary in planes which go through zero, works for nd-arrays
%           'all'                                   Removes all zeros (totally) and returns a vector, works for nd-arrays
%           'left' | 'right' | 'top' | 'bottom'     Works only for 2-arrays. Removes zeros at the given side of the array
%
% Options:
%           'keepdim'                               Does not fill array up with zeros if dimension afterwards is smaller than before
%
% Output:    cout                                   Sequence c, zeros removed
%            origin                                 dim x 1-vector. Idx of first entry of cout, with respect to c. I.e.: cout(idx)=c(1).
%                                                   For 'all' is the linear index of the nonempty entries
%                                                   For 'compact' origin is empty
%                                                   If cout==[], then origin = zeros(dim,1)
%
% E.g.: removezero(padarraym([0 0 0 0;0 1 0 0;0 0 0 1;0 0 0 1],[0 0 1 1]),[1 2 3  ])
%       removezero([0 0 0 0;0 1 0 0;0 0 0 1;0 0 0 1],'border')
%       removezero([0 ],2)
%
% Written by: tommsch, 2018

% Changelog: 2020-04-28, tommsch:       Behaviour change: Row vectors are now considered to be 2dimensional arrays

%XX seems not work for complex entries



c = varargin{1};
dim = ndimsm( c ); %dimension of original array

if( parsem('all',varargin) );
    c = c(:); 
    idx = any( c, 2 );
    origin = find( idx, 1 );
    c( ~idx, : ) = [];
elseif( parsem( 'border',varargin) || parsem( 'outside',varargin) || parsem( 'inside',varargin) );
    [c,origin]=removezeros_border(varargin{:});
elseif( size(varargin,2)>=2 && isnumeric(varargin{2}) ) %removes zeros in directions varargin{2}
    origin = ones( dim, 1 );
    for i = varargin{2}; 
        [c,origin(i)] = removezeros_direction( c, i ); end;
    if( ~all(origin) ); 
        origin = zeros( dim, 1 ); end;
elseif( parsem( 'left',varargin) || parsem( 'right',varargin) || parsem( 'top',varargin) || parsem( 'bottom',varargin) );
    [c,origin] = removezeros_side( varargin{:} );
elseif( parsem( 'compact',varargin) );
    [c,origin] = removezeros_compact( varargin{:} );
else
    error('''removezeros'': wrong argument.'); end;

if(parsem( 'keepdim',varargin) && ~parsem( 'all',varargin));    
    c=makearraygreatagain(c,dim); end;

cout=c;

end

function c = makearraygreatagain( c, targetdim )
currentdim = ndimsm( c ); 
if( currentdim<targetdim )
    val = ones( 1, targetdim );
    val = num2cell( val );
    val{end} = 2;
    c(val{:}) = 0;
    if( targetdim==1 && ~iscolumn(c) ); 
        c=c.'; end; end;
end

function [ c, origin ] = removezeros_direction( c, dir )
    dim = ndimsm( c ); % get dimension
    o = [ 1:dir-1 dir+1:dim]; %all other dimensions
    IDX = cell( 1, dim );
    
    for i = o; 
        IDX{i} = 1:size( c, i ); end;
    IDX{dir} = anym( c, o ); %find zeros
    if( isequal(IDX{dir},0) ); 
        c = []; 
    else; 
        c = c(IDX{:}); end;
    origin = find( IDX{dir} ,1 );
    if( isempty(origin) ); 
        origin = 0; end;
end

function [ c, origin ] = removezeros_border( varargin )
    c = varargin{1};
    dim = ndimsm(c); % get dimension
    %if( isvector(c)==1 );
    %    dim = 1; end;

    inside = parsem( 'inside', varargin );
    outside = parsem( 'outside', varargin );

    IDX = cell( 1, dim );
    origin = zeros( dim, 1 );
    for i = 1:dim;
        o = [ 1:i-1 i+1:dim];
        idx = anym( c, o ); %take slices 
        idx = idx(:);
        FIRST = find( idx, 1, 'first' );
        LAST = find( idx, 1, 'last' );
        if( inside );  
            LAST  = size( c, i ); end;
        if( outside );  
            FIRST = 1; end;
        if( isempty(FIRST) ); 
            origin = 0; 
        else; 
            origin(i) = FIRST; end;
        IDX{i} = FIRST:LAST; end;
    c = c(IDX{:});
end


function [ c, origin ] = removezeros_compact( varargin )
    c = varargin{1};
    dir = parsem( 'compact', varargin, 0 );
    
    if( length(dir)>1 ); 
        for i = 1:length( dir );
            c = removezero( c, 'compact', dir(i) ); end; 
    else
        len = size( c, dir );
        c = num2cell( c, dir );
        nd = ndims( c );
        for i = 1:numel( c );
            c{i} = removezero( c{i}, 'all' ); %remove all zeros
            c{i}(len+1,1) = 0; %add one zero at the end, do not test for length here
            if( dir~=1 ); 
                ordering = 1:nd;
                ordering(1) = dir; 
                ordering(dir) = 1;
                c{i} = permute( c{i}, ordering ); end; end;
        c = cell2mat( c ); 
        c = removezero( c, dir ); end;
    origin = [];
end

function [c,origin]=removezeros_side(varargin)

    left = parsem( 'left', varargin );
    right = parsem( 'right', varargin );
    top = parsem( 'top', varargin );
    bottom = parsem( 'bottom', varargin );


    c = varargin{1};
    dim = ndimsm(c);
    if( ~ismatrix(c) ); 
        error( 'removezero:matrix', 'For this option, array must be a 2-array.' ); end;

        if( left )
            chelp = sum( c~=0, 1 );
            idx = find( chelp, 1 );
            c = c(:,idx:end);
            origin = ones( dim, 1 ); 
            origin(2) = idx;
        elseif( right )
            chelp = sum( c~=0, 1 );
            idx = find( chelp, 1, 'last' );
            c = c(:,1:idx);
            origin = ones( dim, 1 );
        elseif( top )
            chelp = sum( c~=0, 2 );
            idx = find( chelp, 1 );
            c = c(idx:end,:);
            origin = ones( dim, 1 ); 
            origin(1)=idx;
        elseif( bottom )
            chelp = sum( c~=0, 2 );
            idx = find( chelp, 1, 'last' );
            c = c(1:idx,:);
            origin = ones( dim, 1 ); end;
end