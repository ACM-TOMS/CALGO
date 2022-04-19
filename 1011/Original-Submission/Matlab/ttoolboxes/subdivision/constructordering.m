function [ oo] = constructordering(varargin)
% oo = constructordering( oo1, [pp1, oo2, pp2, ... ])
% oo = constructordering( oo )
% Constructs an ordering.
% An ordering is a representation of an infinite periodic sequence.
% It is stored as an 1x2 - cell array, where the first cell is the non-periodic part, and the second cell is the periodic part
% Each cell can have an arbitrary number of rows 
% E.g.: {[2 3],[5 0]} correspdonds to 23505050505050505050
%
% Input:
%   oo_i, pp_i      vector of numbers which are the non-periodic/periodic part for the i-th row
%                   Number of arguments is either 1 or an even number
%   oo              ordering, in this case oo is returned without change
%
% Options:
%   'random',val    Generates a random ordering of values in oo1/pp1
%  
% Output:
%   oo          The ordering defined by the input arguments. 
%
% Note: If oo1 is an ordering, oo1 is returned
%       All periodic and non-periodic parts must have the same length.
%
% Eg: vdisp(constructordering([1 2],[3 5 6],[10 20],[30 50 60]))
%     vdisp(constructordering(2,2,'random',6))
%     
% See also: isordering, vector2ordering, ordering2vector, findperiod
%
% Written by: tommsch, 2018

% XX Change that periodic parts must have the same length

 %#ok<*ALIGN>

if( isordering(varargin{1}) ); 
    if( nargin>1 ); 
        error('An ordering is given, but there are more arguments. This is wrong.'); end;
    oo = varargin{1}; 
    return; end;

[randomlength,varargin] = parsem( 'random', varargin, 0 );


oo = cell(1,2);
if( ~iswholenumber(nargin/2) ); 
    varargin{end+1} = []; end;
NROWS = size( varargin, 2 )/2;

L1 = max( length(varargin{1}), randomlength );
L2 = max( length(varargin{2}), randomlength );
oo{1} = zeros(NROWS,L1);
oo{2} = zeros(NROWS,L2);

for i = 1:NROWS
    val1 = varargin{2*i-1}; 
    val2 = varargin{2*i};
    %if( val1==0 );  %XX Commented out by tommsch 2019_30_10, because I have no idea whats that commands for.
    %    oo{1} = []; end;
    %if( val2==0 ); 
    %    oo{2} = []; end;
    if( randomlength );
        if( val1 ); 
            oo{1}(i,:) = randi( val1, 1, randomlength ); end;
        if( val2 ); 
            oo{2}(i,:) = randi( val2, 1, randomlength ); end;
    else
        oo{1}(i,:) = varargin{2*i-1};
        oo{2}(i,:) = varargin{2*i}; end; end

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 