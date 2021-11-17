function C = times(A,B)
% times(A,B)
% Scalar * Cell / Cell * Scalar: times elementwise
% Matrix * Cell / Cell * Matrix: mtimes elementwise
% Cell * Cell   : times elementwise
%
% Written by: tommsch, 2018

if( iscell(A) && ~iscell(B) && isscalar(B) )
    C = cellfun( @(x) x.*B, A, 'UniformOutput',false );  %# Apply mtimes cell-wise
elseif( iscell(A) && ~iscell(B) )
    C = cellfun( @(x) x*B, A, 'UniformOutput',false );  %# Apply mtimes cell-wise
elseif( ~iscell(A) && isscalar(A) && iscell(B) )    
    C = cellfun( @(x) A.*x, B, 'UniformOutput',false );  %# Apply mtimes cell-wise
elseif( ~iscell(A) && iscell(B) )
    C = cellfun( @(x) A*x, B, 'UniformOutput',false );  %# Apply mtimes cell-wise
elseif( iscell(A) && iscell(B) )
    C = cellfun( @(x,y) x.*y, A, B, 'UniformOutput',false ); end;

end