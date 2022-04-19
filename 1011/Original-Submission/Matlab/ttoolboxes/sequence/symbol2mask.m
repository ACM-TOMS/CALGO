function [ a, amin ] = symbol2mask(varargin)
% [ a, amin ]     = symbol2mask( symbol, [options])
% [ {a}, {amin} ] = symbol2mask( symbol || {symbol}, [options])
% Computes the masks from given symbols.
%
% Input: 
%       symbol              symbol
%       {symbol}            cell array of symbols
%
% Options:
%       'dim',val           specifies the dimension
%       'var',val           Controls the name of the variables. Default: "symvar(symbol{:})"
%                           Format string:  char (ex: 'x')                                  Variables: char,char+1,char+2, (ex: x,y,z)
%                                   char+num  (ex:'b3')                             Variables: charnum,charnum+1,charnum+2 (ex: b3,b4,b5,b6,...)
%                                   cellarray of strings (ex:{'x','a23','fa')      Variables: cellarray{1},cellarray{2},... (ex: x,a23,fa)
%                                   array of syms(ex:[x,a23,f]                      Variables: array(1), array(2), ... (ex: x,a23,fa)
%                                         If one puts values in the cell array, then the symbol is evaluated at this point.
%       'sym'               output is symbolic (default is double)
%       'verbose',val       Verbose level
%
% Output:
%       a                   mask (if "symbol" is given) or
%       amin                corresponding value of amin (see getS)
%   or
%       {a},{amin} as above, cell array with the same size as {symbol}
%
% Note: This function is much faster than mask2symbol
%                  
%
% E.g.: symbol2mask(mask2symbol([1 2 3; 4 5 6]))
%       syms  x y z; symbol2mask(2*z+3*x*z^2,'var',{x y z })
%       syms  x y z; vdisp(symbol2mask({2*x+y,3*x*y^2+4}))
%
% See also: mask2symbol, coeffs, getS
%
% Written by: tommsch, 2017

% XX Function could be rewritten, since it is very ugly

%#ok<*ALIGN>

dim = parsem( 'dim', varargin, 0 );
symflag = parsem( 'sym', varargin );
var = parsem( 'var', varargin, [] );
verbose = parsem( {'verbose','v'}, varargin, 1 );

symbolgiven = 0;
if( issym(varargin{1}) ) %sy is a cell array of masks afterwards
    symbolgiven = 1;
    sy = {varargin(1)};
else
    sy = varargin{1}; end; 
amin = sy; %initialize amin
a = sy;
J = numel( sy ); %number of symbols

%get dimension (and create variables)
sym SYMS;
SYMS = 0;
for j = 1:J
    for k = 1:numel( sy{j} );
        SYMS = sy{j}{k}+SYMS; end; end;

zin = symvar(SYMS);
if( ~dim )
    dim = length( zin ); end;

idx = 1;
z = sym( zeros(1,dim) ); 
if( iscell(var) )
    for k = 1:length( var );
        z(idx) = sym( sprintf('%c', var{k}) ); 
        idx = idx+1; end; 
    dim = length( z ); %set dimension to new value
elseif( issym(var) )
    for k = 1:length( var );
        z(idx) = var(k);
        idx = idx+1; end;
elseif( isempty(var) )
        z = zin;      %keep variables in z
elseif( length(var)==1 );
    for k = double( var(1) ):dim+var(1)-1; 
        z(idx) = sym( sprintf('%c', k) ); 
        idx = idx+1; end; 
elseif( length(var)==2 );
   z = sym( zeros(1,dim) ); 
   for k = str2double( var(2) ):dim+str2double( var(2) )-1; 
       z(idx) = sym( sprintf([num2str(var(1)) '%d'],k) ); 
       idx = idx+1; end; end;

if( dim~=numel(z) ); 
    vprintf( 'Dimension and number of variables are not the same. Output could be wrong. \nGive <''var'',val>  instead of <''dim'',val> to specify the symbols directly instead.\n', 'imp',[1 verbose] ); end;

for j = 1:J;
    for k = 1:numel( sy{j} );
        if( isempty(sy{j}{k}) || isequal(sy{j}{k},sym(0)) || isequal(sy{j}{k},0) );  %if the symbol is empty
            a{j}{k} = 0; 
            amin{j}{k} = zeros( dim, 1 );
        else
            [~,den] = numden( sy{j}{k} );
            denval = coeffs( den );
            mincoeff = coeffs( den, z, 'All' );
            mincoeff = -(size( mincoeff )-1)';
            if( dim==1 ); 
                mincoeff = mincoeff(2); end;
            amin{j}{k} = mincoeff;
            a{j}{k} = coeffs( expand(sy{j}{k}*den/denval) , z, 'All' );
            if( dim==1 ); 
                a{j}{k} = a{j}{k}.'; end; end;
        for l = 1:max( dim, 2 ); 
            a{j}{k} = flip( a{j}{k}, l ); end; %coeffs is ordering the coefficient from high to low, but we are ordering them low to high
        if( ~symflag );  
            a{j}{k} = double( a{j}{k} ); 
            amin{j}{k} = double( amin{j}{k} ); end; end; end;


if( symbolgiven )
    a = a{1}{1};
    amin = amin{1}{1};
else
    a = a{1};
    amin = amin{1}; end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 