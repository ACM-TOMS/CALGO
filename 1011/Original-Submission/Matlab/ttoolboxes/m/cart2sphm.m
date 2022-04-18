function [ sco, typeout ] = cart2sphm( X, type )
% [ sco, typeout ] = cart2sphm( X, [type] )
% Transforms cartesian to hypershperical coordinates
%
% Input:
% =====
%   X       dim X N vector, cartesian coordinates
%
% Options:
% =======
%   type    int, default=0, type of coordinates to transform to (only affects the first two angles)
%           last entry in vector is always the radius \in [0, inf]
%           Positive types are in radians, negative types are in degree
%           some types have multiple identifiers
%
%       Hyperspherical coordinates:
%       ===========================
%           0,1/-1 ('hsph','h'/'hsphdeg','hd')              Hyperspherical coordinates (phi_1, ..., phi_{dim-1}, r)
%                       \phi_1, \phi_2, ... \phi_{dim-2} \in [0, pi]/[0, 180],
%                       \phi_{dim-1} \in [-pi, pi)/[-180, 180),
%                       does not cooincide in dim==3 with matlabs cart2sph
%                       cooincides in dim==2 with polarcoordinates where \phi_1 \in [-pi, pi)/[-180, 180) and
%                                                 matlabs cart2pol, type +-2, type +-3 (up to round-off errors)
%                       In dim==1 these coordinates do not make sense, the function returns the input value, This bevahiour is subject to be changed
%                       (old option 0/10)
%           2/-2   ('azel','ae'/'azeldeg','aed')            Azimuth/Elevation coordinates, only for dim==2,3
%                       \phi_1 \in [-pi,    pi)/[-180, 180), from x-axis to y-axis counter-clockwise, (azimuth for S^2)
%                       \phi_2 \in [-pi/2, pi/2]/[-90,  90], from xy-plane to z-axis, (elevation for S^2) (only for dim==3)
%                       uses matlabs cart2sph
%                       (old option 1/11)
%           3/-3   ('azinc','ai'/'azincdeg','aid')          Azimuth/Inclination coordinates, only for dim==2,3
%                       \phi_1 \in [-pi,  pi)/[-180, 180), from x-axis to y-axis counter-clockwise (azimuth for S^2)
%                       \phi_2 \in   [0,  pi]/[   0, 180], from z-axis to xy-plane (only for dim==3)
%                       uses matlabs cart2sph
%                       (old option 2/12)
%           4/-4   ('ant'/'antdeg','antd')                  Antenna pattern coordinates, spherical coordinates on the xy-plane and the yz-plane, only for dim==2,3
%                       \phi_1 \in [0, 2*pi]/[0, 360], from y-axis (North) to x-axis clockwise
%                       \phi_2 \in [0, 2*pi]/[0, 360], from y-axis (North) to minus-z-axis clockwise 
%                       only allowed in sph2cartm, sph2cartm, sph2sphm, sph2sphm2 as type1, 
%                       (old option -1/9)
%
%       Hypercylindrical coordinates:
%       =============================
%           100,101/-100/-101 ('hcyl','hc'/'hcyldeg','hcd')     Hypercylindrical coordinates in ordering (phi, x_3, x_4, ..., x_{dim}, r)
%                                                               Cartesian product of cylinder with line segments
%                       \phi \in [-pi, pi)/[-180, 180), from x_1-axis (x-axis) to x_2-axis (y-axis) counter-clockwise
%                       Note the strange ordering of coordinates!
%                       In dim==3 these coordinate cooincide with Spherinder coordinates
%                       In dim==2 Hypercylindrical and Hyperspherical (i.e. polar) coordinates cooincide
%                       In dim==1 these coordinates do not make sense
%           104/-104 ('hant'/'hantdeg','hantd')                 Strange cylindrical coordinates in ordering (phi,z,r), only for dim==3, 
%                       \phi \in [0, 2*pi]/[0, 360], from y-axis (North) to x-axis clockwise
%                       only allowed in sph2cartm, sph2cartm, sph2sphm, sph2sphm2 as type1
%                       Note the strange ordering of coordinates!
%
%      Not implemented yet:
%      ====================
%           105/-105 ('spherinder'/'spherinderdeg','spherinderd')           Spherinder coordinates in ordering (phi_1, phi_2, ..., phi_{dim-2}, x_{dim}, r)
%                                                                           Cartesian product of hypersphere with line segment
%                       \phi_1, \phi_2, ... \phi_{dim-3} \in [0, pi]/[0, 180],
%                       \phi_{dim-2} \in [-pi, pi)/[-180, 180),
%                       \phi \in [0, 2*pi)/[0, 360), from x_1-axis (x-axis) x_2-axis (y-axis) counter-clockwise
%                       Note the strange ordering of coordinates!
%                       In dim==3 these coordinate cooincide with Hypercylindrical coordinates
%                       In dim==2 these coordinates do not make sense
%                       In dim==1 these coordinates cooincide with Cartesian coordinates
%
%
% Experimental Options:
% =====================
%   The string aliases are experimental
%   The short string aliases are subject to be changed in future versions
%   Use the short versions with causions, since they may collide in some functions with other short versions of commands
%
% Output:
% =======
%   sco         angles of spherical coordinates, last entry: radius
%   typeout     the deduced type
%
% Note:     This function, and the corresponding ones given below, partially work for dim==0 and dim==1.
%           Some types only work for certain dimensions
%           Due to roundoff errors, the angles may lie slightly outside their documented range.
%
% E.g.: cart2sphm( [2 3 1]' )
%
%
% See also: sph2cartm, cart2sphm2, sph2cartm2, sph2sphm, sph2sphm2
%
% Written by: tommsch, 2019


% Changelog: tommsch, 2019-11-25,   sco can be a matrix of column-vectors, also in all corresponding functions
%            tommsch, 2019-11-29,   Added ''type''
%            tommsch, 2020-03-06,   Added option 20/-20, 
%                                   Changed values of types
%                                   Added experimental aliases for types
%                                   Added return of type
%                                   Added experimental support of type 2,3 for dim==2



% XX Implement types
%                   spherical, only possible for dim>=3
%                   \phi_1 \in [0,2\pi), from x-axis to y-axis counter-clockwise, (azimuth for S^2)
%                   \phi_{dim-1} \in [-\pi,\pi], (elevation for S^2)
%                   cooincides in dim=3 with matlabs cart2sph
% XX Implement spherinder
% XX Implement hyperspherical for any dimensions
% XX Implement Longitutde/Latitude coordinates

%parse input
if( nargin<=1 || isempty(type) );
    type = 1; end;

if( ischar(type) )
    switch type
        case {'hsph','h',''};       type =    0;
        case {'hsphdeg','hd'};      type =   -1;
        case {'azel','ae'};         type =    2;
        case {'azeldeg','aed'};     type =   -2;
        case {'azinc','ai'};        type =    3;
        case {'azincdeg','aid'};    type =   -3;
        case {'ant'};               type =    4;
        case {'antdeg','antd'};     type =   -4;
        case {'hcyl','hc'};         type =  100;
        case {'hcyldeg','hcd'};     type = -100;
        case {'hant'};              type =  104;
        case {'hantdeg','hantd'};   type = -104;
        otherwise; error( 'Wrong string for type.' ); end; end;
    
typeout = type;

if( type==9 || type==10 || type==11 || type==12 )
    fprintf('Old type used. See the help of this function for the new numbers of the types.\nFix calling arguments.\n'); error('Error in cart2sphm'); end;

% handle abnormal cases
if( nargin==0 || isempty(X) ); 
    sco=[]; 
    return; end;

dim = size( X, 1 );
if( dim==1 ); 
    sco=X;
    return; end;

%transform coordinates and do checks
switch type
    case {0,1,-1}; %old case 0/10
        sco = zeros( dim, size(X,2) );
        for i = 1:dim-1;
            sco(i,:) = acos( X(i,:)./sum(X(i:end,:).^2,1).^(1/2) ); end;
        idx = X(end,:)<0;
        sco(dim-1,idx) = -sco(dim-1,idx);
        sco(dim,:) = sum(X.^2,1).^(1/2);
        idx = isnan(sco);
        sco(idx) = 0;
        sco(isnan(sco)) = 0; 
        
    case {2,-2,3,-3};  %old case 1,2,11,12
        switch dim
            case 2;
                [th,rho] = cart2pol( X(1,:), X(2,:) ); %transforms to [-pi,pi]
                sco = [th;rho];
            case 3;
                [az,el,r] = cart2sph( X(1,:), X(2,:), X(3,:) );
                sco = [az;el;r];
            otherwise
                error( '''type +-2/3'' only implemented for dim==2,3 yet.' ); end;
        
    case {100,101,-100,-101};
        if(dim==1);
            error( '''Type +-100, +-101'' (Hyperspherical coordinates) do not make sense for dim==1.' );  end;
        [th,rho] = cart2pol( X(1,:), X(2,:) ); 
        sco = [th; X(3:end,:); rho];
        
    case {4,104,-4,-104};
        error( 'Not possible to transform to ''type +-4, +-104'' (Antenna pattern) coordinates.' )
        
    otherwise;
        error( 'Unkown type.' ); end;
    
%post-processing    
switch type
    case {3,-3}; %old case 2,12
        switch dim
            case 2;
                %do nothing
            case 3;
                sco(2,:) = pi/2-sco(2,:);
            otherwise;
                error('Fatal error. Check Code.'); end;
    otherwise; 
        %fallthrough
    end; 
        
%transform values to degrees if necessary
switch type
    case {0,1,2,3,4,100}; %old case 0,1,2,-1
        %do nothing
    case {-1,-2,-3,-4}; %old case 9,10,11,12
        sco(1:end-1,:) = sco(1:end-1,:)*180/pi;
    case {-100,-101};
        sco(1,:) = sco(1,:)*180/pi;
    case {-105};
        sco(1:end-2,:) = sco(1:end-2,:)*180/pi;
    otherwise;
        error( 'Unkown type.' ); end;    

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   