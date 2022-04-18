function [ X, typeout ] = sph2cartm2(varargin)
% [ X, typeout ] = sph2cartm2( sco, [type] )
% Transforms hypershperical to cartesian coordinates, assumes that the radius is one.
% i.e.: size(sco)==[(dim - 1), N] yields size(X)==[dim,N]
%
% See also: cart2sphm
%
% Written by: tommsch, 2019

sco = varargin{1};
[X,typeout] = sph2cartm( [sco; ones(1,size(sco,2))], varargin{2:end} );

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   