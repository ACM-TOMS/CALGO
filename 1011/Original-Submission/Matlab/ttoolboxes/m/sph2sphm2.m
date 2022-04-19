function [ sco, typeout1, typeout2 ] = sph2sphm2(sco, type1, type2, flag)
% [ sco, typeout1, typeout2 ] = sph2sphm2( sco, type1, type2, [flag] )
% [ sco, typeout1, typeout2 ] = sph2sphm2( sco, type )
% Transforms hypershperical coordinates (of type1) to hyperspherical coordinates (of type2), assumes that the radius is one.
% i.e. size(sco)=[dim-1, N]
% Input:
%   sco         spherical coordinates, 
%                   default: type1 = type2 = 0 (if type1 and type2 are not given or empty)
%                            type2 = type1     (if type2 is not given or empty)
%   flag       default=0
%              if given, disregards the last row in sco while computing and appends it afterwards again.
%              if given, type1 and type2 must be given too
%
% Note: This function supports type "+-4" as type1
%       This function may change the values of sco even if type1==type2
%
% E.g.: 
%
% See also: cart2sphm
%
% Written by: tommsch, 2018

if( isempty(sco) );
    return; end;

%parse input
if( nargin<=1 || (nargin>=2 && isempty(type1)) ); 
    type1 = 0; 
    type2 = 0; end;
if( nargin<=2 || (nargin>=3 && isempty(type2)) );
    type2 = type1; end;
if( nargin<=3 || nargin>=4 && isempty(flag) );
    flag = 0; end


%transform coordinates
if( flag );
    lastrow = sco(end,:); 
    sco=sco(1:end-1,:); end;

if( type1==-3 && type2==-3 && size(sco,2)==2 ) 
    %fast transformation of type -3 (az/inc in degree) in dimension 3
    %phi
    sco(2,:) = mod( sco(2,:), 360 );
    idx = sco(2,:)>=180;
    sco(1,idx) = sco(1,idx)+180;
    sco(2,idx) = 360-sco(2,idx);
    %theta
    sco(1,:) = mod( sco(1,:), 360 );
    idx = sco(1,:)>=180;
    sco(1,idx) = sco(1,idx)-360;
    %typeout
    typeout1=type1;
    typeout2=type2;
else
    [X,typeout1] = sph2cartm2( sco, type1 );
    [sco,typeout2] = cart2sphm2( X, type2 ); end;
    
if( flag );
    sco = [sco; lastrow]; end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   