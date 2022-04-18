%% pntsGeometry
% Coordinates of points for given geometries in two dimensional space.
%
%% Syntax
%
%   [Pnts,Nb,dS] = pntsGeometry(seti,type,varargin)
%
%   [Pnts,Nb,dS] = pntsGeometry(seti,'manually',Pnts)
%   [Pnts,Nb,dS] = pntsGeometry(seti,'circle',Nb,rad)
%   [Pnts,Nb,dS] = pntsGeometry(seti,'square',NbEdge,edgeLength)
%   [Pnts,Nb,dS] = pntsGeometry(seti,'line',NbLine,lineLength,linePos)
%   [Pnts,Nb,dS] = pntsGeometry(seti,'borehole',boreNbLine,boreLineLength,boreLinePos,incMeas)
%
%% Description
% |[Pnts,Nb,dS] = pntsGeometry(seti,type,varargin)|
% computes the coordinates in |Pnts| and serves their number |Nb| and 
% an approximation of the infinitesimal element of closed contour in |dS|.
% The geometry type is defined in |type| and the details in |varargin|
% depending on the type. |seti| is only used to check that |seti.dim = 2|.
%
% |[Pnts,Nb,dS] = pntsGeometry(seti,'manually',Pnts)|
% counts the number |Nb| of input points |Pnts| and sets |dS| to 1. Input
% and output |Pnts| is the same.
%
% |[Pnts,Nb,dS] = pntsGeometry(seti,'circle',Nb,rad)|
% generates coordinates |Pnts| on a cricle. The number of |Nb| points are 
% equidistant arranged on a circle with radius |rad|. 
% Distance |dS| is computed.
%
% |[Pnts,Nb,dS] = pntsGeometry(seti,'square',NbEdge,edgeLength)|
% generates coordinates |Pnts| on a square. The number of |NbEdge|+1 points
% are equidistant arranged on one edge with length |edgeLength|. 
% (Total number of points in square is |4*NbEdge|.
%  |NbEdge|+1 per edge is because the points in corners belonging to two corners.)
% This is done for all four edges of the square. Distances |dS| are computed.
%
% |[Pnts,Nb,dS] = pntsGeometry(seti,'line',NbLine,lineLength,linePos)|
% generates |NbLine| coordinates |Pnts| on a line with length |lineLength|
% with a position of |linePos| as distance to origin. |dS| is set to 1.
%
% |[Pnts,Nb,dS] = pntsGeometry(seti,'borehole',boreNbLine,boreLineLength,boreLinePos,incMeas)|
% generates |boreNbLine| coordinates |Pnts| on a line with length |boreLineLength|
% with a position of |boreLinePos| as distance to origin. |dS| is set to 1.
% Set |incMeas| to |'inc'| or |'meas'|, then |Pnts| are shifted such that
% points of both types will be equidistant.
%
%% Examples
%
% *Example 1: manually*
%
% 3 manually positioned points with coordinates (1,0), (0.5,1.5), and (-0.5,1.5).
%
%   seti.dim = 2;
%   Pnts = [1 0.5 -0.5;
%           0 1.5  1.5];
%   [Pnts,Nb,dS] = pntsGeometry(seti,'manually',Pnts);
%   figure(101); plot(Pnts(1,:),Pnts(2,:),'b.','MarkerSize',30); axis square;
%
% <<../extGraph/pntsGeometry_fig101.png>>
%
%
% *Example 2: circle*
%
% Arrange 10 points on a circle with radius 5.
%
%   seti.dim = 2;
%   Nb = 10;
%   rad = 5;
%   [Pnts,Nb,dS] = pntsGeometry(seti,'circle',Nb,rad);
%   figure(102); plot(Pnts(1,:),Pnts(2,:),'b.','MarkerSize',30); axis square;
%
% <<../extGraph/pntsGeometry_fig102.png>>
%
%
% *Example 3: square*
%
% Arrange 4+1 points on each edge of a square with length 3.
% Total number of points is 4*4 = 16.
%
%   seti.dim = 2;
%   NbEdge = 4;
%   edgeLength = 3;
%   [Pnts,Nb,dS] = pntsGeometry(seti,'square',NbEdge,edgeLength);
%   figure(103); plot(Pnts(1,:),Pnts(2,:),'b.','MarkerSize',30); axis square;
%
% <<../extGraph/pntsGeometry_fig103.png>>
%
%
% *Example 4: line*
%
% Arrange 5 points on a line with length 3, right to origin with
% distance 2.
%
%   seti.dim = 2;
%   NbLine = 5;
%   lineLength = 3;
%   linePos = 2;
%   [Pnts,Nb,dS] = pntsGeometry(seti,'line',NbLine,lineLength,linePos);
%   figure(104); plot(Pnts(1,:),Pnts(2,:),'b.','MarkerSize',30); axis square;
%
% <<../extGraph/pntsGeometry_fig104.png>>
%
%
% *Example 5: borehole*
%
% Arrange 10 transmitters and 10 receivers on a line of length 4 with a
% position 5 right from origin.
%
%   seti.dim = 2;
%   boreNbLine = 10;
%   boreLineLength = 4;
%   boreLinePos = 5;
%   [incPnts,incNb,dSInc] = pntsGeometry(seti,'borehole',boreNbLine,boreLineLength,boreLinePos,'inc');
%   [measPnts,measNb,dSMeas] = pntsGeometry(seti,'borehole',boreNbLine,boreLineLength,boreLinePos,'meas');
%   figure(105); 
%   hold on;
%   plot(incPnts(1,:),incPnts(2,:),'b.','MarkerSize',30);
%   plot(measPnts(1,:),measPnts(2,:),'rs','MarkerSize',7,'MarkerFaceColor','red');
%   hold off;
%   axis square;
%
% <<../extGraph/pntsGeometry_fig105.png>>
%
%
%% Input Arguments
%
% * seti        :   structure array which has to contain |seti.dim = 2|
%                   because two dimensional space
% * type        :   string with type of geometry:
%                   'manually', 'circle', 'square', 'line', 'borehole'.
% * varargin    :   various input arguments depending on choosen type.
%
% *varargin in case of type = 'manually'*
%
% * Pnts        :   coordinates of points (real matrix of size 2 x Nb).
%
% *varargin in case of type = 'circle'*
%
% * Nb          :   Number of points
% * rad         :   radius of circle
%
% *varargin in case of type = 'square'*
%
% * NbEdge      :   NbEdge+1 points are on one edge
%                   (+1 because a point in a corner belongs to two corners).
%                   (Total number of points in square is |4*NbEdge|.)
% * edgeLength  :   Length of one edge.
%
% *varargin in case of type = 'line'*
%
% * NbLine      :   Number of points on line.
% * lineLength  :   Length of line.
% * linePos     :   Distance of line right to origin.
%
% *varargin in case of type = 'borehole'*
%
% * boreNbLine      :   Number of points on line.
% * boreLineLength  :   Length of line.
% * boreLinePos     :   Distance of line right to origin.
% * incMeas         :   Type of points, i.e. 'inc' or 'meas', to shift the
%                       points, see Example 5 above.
%
%
%% Output Arguments
%
% * Pnts    :   Coordinates of points (real matrix of size 2 x Nb).
% * Nb      :   Number of points.
% * dS      :   Approximation of the infinitesimal element of closed 
%               contour with control points.
%               A closed contour does not make sense in the cases
%               'manually', 'line', and 'borehole', then it is set to 1.
%
%
%% More About
%
% * currently only in two dimensional space
%   (this is the reason for seti.dim = 2)
%
%
%% See Also
%
% * <dS2D.html>
% * <expSetup.html>
% * <setIncPnts.html>
% * <setMeasPnts.html>
%
%
%% Code: pntsGeometry
%
function [Pnts,Nb,dS] = pntsGeometry(seti,type,varargin)

if seti.dim ~= 2
    error('pntsGeometry currently needs 2D.')
end

%%
% *Split varargin in dependence of parameter type*

if strcmp(type,'manually')
    %   pntsGeometry(seti,'manually',Pnts)
    Pnts = varargin{1};
    
elseif strcmp(type,'circle')
    %   pntsGeometry(seti,'circle',Nb,rad)
    Nb = varargin{1};
    rad = varargin{2};
    Pnts = circle(Nb,rad);
    
elseif strcmp(type,'square')
    NbEdge = varargin{1}; % number of points on one edge
    edgeLength = varargin{2}; % length of edge
    %   pntsGeometry(seti,'square',NbEdge,edgeLength)
    Pnts = square(NbEdge,edgeLength);
    
elseif strcmp(type,'line')
    NbLine = varargin{1};
    lineLength = varargin{2};
    linePos = varargin{3}; % position of line (distance to origin)
    %   pntsGeometry(seti,'line',NbLine,lineLength,linePos)
    Pnts = line(NbLine,lineLength,linePos);
    
elseif strcmp(type,'borehole')
    NbLine = varargin{1};
    lineLength = varargin{2};
    linePos = varargin{3}; % position of line (distance to origin)
    incMeas = varargin{4}; % 'inc' or 'meas'
    %   pntsGeometry(seti,'borehole',NbLine,lineLength,linePos,incMeas)
    Pnts = borehole(NbLine,lineLength,linePos,incMeas);
end

%% 
% *Compute Nb and dS*

[~,Nb] = size(Pnts); % in some cases (circle, line) it is input, but not in all...
% Note: Nb = length(Pnts); This is not correct in case of seti.incNb = 1 it
% results in 2, which is wrong...; that is why size is used...

if strcmp(type,'manually') || strcmp(type,'line') || strcmp(type,'borehole')
    dS = 1;
else
    % transpose of Pnts because dS2D expects Nb x 2 and not 2 x Nb
    dS = dS2D(transpose(Pnts));
end

end

%% Code: Subfunctions

%%
%
% *circle*
%
function Pnts = circle(Nb,rad)
% set up uniformly spaced points
%
% * Nb  :   number of points on circle
% * rad :   radius
%
cmplxPnts = exp(2*pi*1i*(0:1:Nb-1)./Nb);
Pnts = [real(rad.*cmplxPnts); imag(rad.*cmplxPnts)];
%     % next line gives example for measurements on circle on the right of the scattering medium
%     seti.incPnts = [real(0.3+seti.radSrc.*cmplxPnts); imag(0.3+seti.radSrc.*cmplxPnts)];
end

%%
% *square*
%
function Pnts = square(NbEdge,edgeLength)

dist = edgeLength/2;
% add 2 and then delete the corner points
x = linspace(-dist,+dist,NbEdge+1);
x = x(1:end-1); % without last entry

distVec = dist*ones(size(x));

% edges
top    = [distVec; x];
bottom = [-distVec; -x];
left   = [x; -distVec];
right  = [-x; +distVec];

Pnts = [top right bottom left];
% Pnts are arranged in mathematical negative sense starting top left.
% Make sure that the vectors are stored in a suitable order because the distances computed in dS2D.
end

%%
% *line*
function Pnts = line(NbLine,lineLength,linePos)

dist = lineLength/2;
x = linspace(-dist,+dist,NbLine);
Pnts = [linePos*ones(size(x)); x];
end

%%
% *borehole*
%
function Pnts = borehole(Nb,lineLength,linePos,type)
% type: 'inc' or 'meas'

d = lineLength/(2*Nb); % distance between transmitters and receivers
dist = lineLength/2;
if strcmp(type,'inc')
    x = linspace(-dist,+dist-d,Nb);
else % type = meas
    x = linspace(-dist+d,+dist,Nb);
end
Pnts = [linePos*ones(size(x)); x];

end
