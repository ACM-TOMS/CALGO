function [ sco, typeout1, typeout2 ] = sph2sphm(sco, type1, type2)
% [ sco, typeout1, typeout2 ] = sph2sphm( sco, [type1], [type2] )
% Transforms hyperspherical coordinates (of type 1) to hyperspherical coordinates (of type 2)
% Input:
%   sco         hyperspherical coordinates, 
%                   default: type1 = type2 = 0;
%                            type2 = type1     
%
% Note: This function may change the values of sco even if type1==type2
%
% See also: cart2sphm
%
% Written by: tommsch, 2019

if( isempty(sco) );
    return; end;

if( nargin==1 ); 
    type1 = 0; 
    type2 = 0; end;
if( nargin==2 );
    type2 = type1; end;

[X,typeout1] = sph2cartm( sco, type1 );
[sco,typeout2] = cart2sphm( X, type2 );

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   