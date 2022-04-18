%% setIncPnts
% Set the positions of the transmitters.
%
%% Syntax
%
%   seti = setIncPnts(seti)
%   seti = setIncPnts(seti,dispDepth)
%
%% Description
% |seti = setIncPnts(seti)| sets the positions |seti.incPnts| of the 
% transmitters in dependence of input in struct seti. The number of
% transmitters is stored in |seti.incNb|. The approximation of the 
% infinitesimal element of closed contour with control points is stored in
% |seti.dSInc|.
%
% |seti = setIncPnts(seti,dispDepth)| does the same but allows to control
% the displayed messages by |dispDepth|.
%
% * If a closed contour does not make sense, |seti.dSInc| is set to 1.
%
% *In case of 2D* (seti.dim == 2)
%
% * If *point sources* are used (seti.incType = 'pointSource'), all
% geometries in seti.incPntsType can be used 
% (i.e. 'manually', 'circle', 'square', 'line', 'borehole').
% * If *plane waves* are used (seti.incType = 'planeWave'), only geometry
% circle is supported (circle with radius 1 is automatically used in case of plane wave).
%
% *In case of 3D* (seti.dim == 3)
%
% * Only a ball is supported in case of point sources or plane waves (in
% latter case, the radius is set to 1).
%
%
%% Example
%
% Arrange 5 point sources on a circle with radius 3. Coordinates are stored
% in seti.incPnts.
%
%   init;
%   seti.dim = 2;
%   seti.incType = 'pointSource';
%   seti.incPntsType = 'circle';
%   seti.incNb = 5;
%   seti.radSrc = 3;
%   seti = setIncPnts(seti);
%
%% Input Arguments
%
% * seti.dim            :   dimension of the problem (2 or 3)
% * seti.incType        :   Type of incident field: 'pointSource' (default) or 'planeWave'
%
% *Input Arguments in case of 2D (seti.dim = 2) and point sources (seti.incType = 'pointSource')*
%
% * seti.incPntsType    :   string with type of geometry:
%                           'manually', 'circle' (default), 'square', 'line', 'borehole'
%
% It follows a *list of type-depending input parameters* to describe the
% details of the geometry. The names are analog to the parameters in 
% <pntsGeometry.html> (mostly with the prefix inc and a big letter,
% exceptions are marked with emphasized).
% Examples are in <pntsGeometry.html>.
% If these parameters are not set, default values are set automatically
% (see inside the code).
%
% _varargin in case of seti.incPntsType = 'manually'_
%
% * seti.incPnts        :   coordinates of points (real matrix of size 2 x Nb).
%
% _varargin in case of seti.incPntsType = 'circle'_
%
% * seti.incNb             :   Number of points (transmitters)
% * _seti.radSrc_          :   radius of circle
%
% _varargin in case of seti.incPntsType = 'square'_
%
% * seti.incNbEdge      :   incNbEdge+1 points are on one edge
%                           (+1 because a point in a corner belongs to two corners).
%                           (Total number of points in square is |4*incNbEdge|.)
% * seti.incEdgeLength  :   Length of one edge.
%
% _varargin in case of seti.incPntsType = 'line'_
%
% * seti.incNbLine      :   Number of points on line.
% * seti.incLineLength  :   Length of line.
% * seti.incLinePos     :   Distance of line right to origin.
%
% _varargin in case of seti.incPntsType = 'borehole'_
%
% * _seti.boreNbLine_      :   Number of points on line.
% * _seti.boreLineLength_  :   Length of line.
% * _seti.boreLinePos_     :   Distance of line right to origin.
%
% *Input Arguments in case of 3D (seti.dim = 3)*
%
% * seti.incPntsType    :   string with type of geometry:
%                           'manually', 'sphereLatLon', 'sphereFibo' (default).
%
% It follows a *list of type-depending input parameters* to describe the
% details of the geometry. The names are analog to the parameters in 
% <pntsGeometry3D.html> (mostly with the prefix inc and a big letter,
% exceptions are marked with emphasized).
% Examples are in <pntsGeometry3D.html>.
% If these parameters are not set, default values are set automatically
% (see inside the code).
%
% _varargin in case of seti.incPntsType = 'manually'_
%
% * seti.incPnts        :   coordinates of points (real matrix of size 3 x Nb).
%
% _varargin in case of seti.incPntsType = 'sphereLatLon' or 'sphereFibo'_
%
% * seti.incNb             :   Number of points (transmitters)
% * _seti.radSrc_          :   radius of sphere
%
% *Optional Input Argument*
%
% * dispDepth   : Depth of displayed messages (0: no, 1 or greater: yes).
%               
%% Output Arguments
%
% * seti.incPnts    :   Coordinates of points (real matrix of size seti.dim x Nb).
% * seti.incNb      :   Number of points.
% * seti.dSInc      :   Approximation of the infinitesimal element of closed 
%                       contour with control points.
%                       A closed contour does not make sense in the cases
%                       'manually', 'line', and 'borehole', then it is set to 1.
%
% For further information see the corresponding parameters 
% Pnts, Nb, and dS in <pntsGeometry.html> and <pntsGeometry3D.html>.
%
%
%% See Also
%
% * <expSetup.html>
% * <pntsGeometry.html>
% * <pntsGeometry3D.html>
% * <dS2D.html>
% * <dS3D.html>
% * <setMeasPnts.html>
%
%
%% Code: function: setIncPnts
%
function seti = setIncPnts(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

if seti.dim == 2
    
    seti = checkfield(seti,'incPntsType','circle',dispDepth); % default
    
    if strcmp(seti.incType,'planeWave')
        % currently only type circle
        [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry(seti,'circle',seti.incNb,1); % radSrc: 1
        
    elseif strcmp(seti.incType,'pointSource')
        switch seti.incPntsType
            case 'manually'
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry(seti,'manually',seti.incPnts);
            case 'circle'
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry(seti,'circle',seti.incNb,seti.radSrc);
            case 'square'
                seti = checkfield(seti,'incNbEdge',9);
                seti = checkfield(seti,'incEdgeLength',7);
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry(seti,'square',seti.incNbEdge,seti.incEdgeLength);
            case 'line'
                seti = checkfield(seti,'incNbLine',35);
                seti = checkfield(seti,'incLineLength',10);
                seti = checkfield(seti,'incLinePos',5);
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry(seti,'line',seti.incNbLine,seti.incLineLength,seti.incLinePos);
            case 'borehole'
                % For type 'borehole' options are the same for inc and meas Pnts.
                % Note the last input argument: inc or meas (to shift)
                seti = checkfield(seti,'boreNbLine',35);
                seti = checkfield(seti,'boreLineLength',500); % 100 seems to low and 1000 seems too high
                % boreLineLength = 500 does not result in a nice reconstruction; change it...
                seti = checkfield(seti,'boreLinePos',5);
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry(seti,'borehole',seti.boreNbLine,seti.boreLineLength,seti.boreLinePos,'inc');
            otherwise
                error('No valid incPntsType.')
        end
        
    else
        error('incType must be planeWave or pointSource.')
    end

elseif seti.dim == 3
    
    seti = checkfield(seti,'incPntsType','sphereFibo',dispDepth); % default
    % seti.incPntsMode = 'incFullSphere'; % old for sphereLatLon
    % sphereLatLon was default until August 2017.
    
    if strcmp(seti.incType,'planeWave')

        switch seti.incPntsType
            case 'sphereLatLon'
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry3D(seti,'sphereLatLon',seti.incNb,1,dispDepth);
            case 'sphereFibo'
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry3D(seti,'sphereFibo',seti.incNb,1,dispDepth);
            otherwise
                error('No valid incPntsType.')
        end
    
    elseif strcmp(seti.incType,'pointSource')

        switch seti.incPntsType
            case 'manually'
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry3D(seti,'manually',dispDepth);
            case 'sphereLatLon'
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry3D(seti,'sphereLatLon',seti.incNb,seti.radSrc,dispDepth);
            case 'sphereFibo'
                [seti.incPnts,seti.incNb,seti.dSInc] = pntsGeometry3D(seti,'sphereFibo',seti.incNb,seti.radSrc,dispDepth);
            otherwise
                error('No valid incPntsType.')
        end
    
    end

else
        disp('   Error: seti.dim is not equal to two or three.')
end

end
