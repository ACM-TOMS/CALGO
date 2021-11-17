function [ sco, typeout ] = cart2sphm2( varargin )
%[ sco, typeout ] = cart2sphm2( X, [type] )
% Transforms cartesian to hyperspherical coordinates not returning the radius.
% i.e.: size(X)==[dim,N] cart2sphm2 yields size(sco)==[(dim - 1), N] 
%
% See also: cart2sphm
%
% Written by: tommsch, 2018

[sco,typeout] = cart2sphm( varargin{:} );
sco(end,:) = [];

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   