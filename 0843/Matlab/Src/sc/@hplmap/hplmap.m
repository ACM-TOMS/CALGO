function map = hplmap(varargin)
%HPLMAP Schwarz-Christoffel half-plane map object.
%   HPLMAP(P) constructs a Schwarz-Christoffel half-plane map object for
%   the polygon P. The parameter problem is solved using default options
%   for the prevertices and the multiplicative constant.
%   
%   HPLMAP(P,OPTIONS) uses an options structure of the type created by
%   SCMAPOPT in solving the parameter problem.
%   
%   HPLMAP(P,Z) creates a hplmap object having the given prevertices Z
%   (the mulitiplicative constant is found automatically).
%   HPLMAP(P,Z,C) also uses the given constant. An OPTIONS argument can
%   be added, although only the error tolerance will be used.
%   
%   HPLMAP(M), where M is a hplmap object, just returns M.
%   
%   HPLMAP(M,P) returns a new hplmap object for the polygon P using the
%   options in hplmap M. The prevertices of M will be used as the
%   starting guess for the parameter problem of the new map. Thus P
%   should properly be a perturbation of the polygon for M. An OPTIONS
%   structure may also be given, to override the options in M.
%   
%   HPLMAP(Z,ALPHA) creates a map using the given prevertices and the
%   interior polygon angles described by ALPHA (see POLYGON help). The
%   image polygon is deduced by computing S-C integrals assuming a
%   multiplicative constant of 1. HPLMAP(Z,ALPHA,C) uses the given
%   constant instead.
%   
%   See also SCMAPOPT, classes POLYGON, SCMAP.

%   Copyright 1998-2001 by Toby Driscoll.
%   $Id: hplmap.m,v 2.5 2002/10/23 18:19:50 driscoll Exp $

superiorto('double');

% For class name with no arguments, return an empty object
if nargin == 0
  map.prevertex = [];
  map.constant = [];
  map.qdata = [];
  map.accuracy = [];
  parent = scmap;
  map = class(map,'hplmap',parent);
  return
end

% Initialize with empties
poly = [];
alpha = [];
z = [];
c = [];
opt = [];
qdata = [];

% Branch based on class of first argument
switch class(varargin{1})
  case 'hplmap'
    map = varargin{1};    
    if nargin == 1
      % Self-return
      return
    else
      % Continuation of given map to given polygon
      poly = varargin{2};
      opt = scmapopt(map);
      z0 = prevertex(map);
      if length(z0) ~= length(poly)
        msg = 'Polygon %s must have the same length as that in %s.';
        error(sprintf(msg,inputname(2),inputname(1)))
      end
      if nargin > 2
        opt = scmapopt(opt,varargin{3});
      end
      opt = scmapopt(opt,'initial',z0);
    end
    
  case 'polygon'
    poly = varargin{1};
    % Parse optional arguments
    for j = 2:length(varargin)
      arg = varargin{j};
      % Each arg is an options struct, z, or c
      if isa(arg,'struct')
        opt = arg;
      elseif length(arg) == length(poly)
        z = arg;
        z = z(:);
      elseif length(arg) == 1
        c = arg;
      else
        msg = 'Unable to parse argument ''%s''.';
        error(sprintf(msg,inputname(j+1)))
      end
    end
    
  case 'double'
    % Args are the prevertex vector, then angle vector
    z = varargin{1}(:);
    alpha = varargin{2}(:);
    if ~isinf(z(end))
      z = [z;Inf];
      alpha = [alpha;1];
    end
    poly = polygon(NaN*alpha*i,alpha);  %  nonsense vertices
    c = 1;
    for j = 3:length(varargin)
      if isa(varargin{j},'struct')
        opt = varargin{j};
      elseif length(varargin{j})==1
        c = varargin{j};
      else
        msg = 'Unable to parse argument ''%s''.';
        error(sprintf(msg,inputname(j+1)))
      end
    end
    
  otherwise
    msg = 'Expected ''%s'' to be a polygon, hplmap, or prevertex vector.';
    error(sprintf(msg,inputname(1)))
    
end % switch


% Retrieve options
opt = scmapopt(opt);

% Take actions based on what needs to be filled in

if isempty(z)
  [w,beta] = scfix('hp',vertex(poly),angle(poly)-1);
  poly = polygon(w,beta+1);

  [z,c,qdata] = hpparam(w,beta,opt.InitialGuess,opt);
end  

if isempty(qdata)
  % Base accuracy of quadrature on given options
  nqpts = ceil(-log10(opt.Tolerance));
  alpha = angle(poly);
  qdata = scqdata(alpha(1:end-1)-1,nqpts);
end

if isempty(c)
  % Find constant
  w = vertex(poly);
  beta = angle(poly)-1;
  idx = 1 + min(find(~isinf(z(2:end))));
  mid = mean(z([1 idx])) + i*diff(real(z([1 idx])))/2;
  I = hpquad(z(1),mid,1,z(1:end-1),beta(1:end-1),qdata) - ...
      hpquad(z(idx),mid,idx,z(1:end-1),beta(1:end-1),qdata);
  c = diff(w([1 idx]))/I;
end

% Get data for the low-level functions
w = vertex(poly);
n = length(w);
beta = angle(poly) - 1;

map.prevertex = z;
map.constant = c;
map.qdata = qdata;

% Make a parent scmap object
parent = scmap(poly,opt);

% Leave a spot for accuracy and create object
map.accuracy = [];
if ~isa(map,'hplmap')
  map = class(map,'hplmap',parent);
else
  map.scmap = parent;
end

% If the polygon was not known, find it from the map
if any(isnan(vertex(poly)))
  poly = forwardpoly(map);
  map.scmap = scmap(poly,opt);
end

% Now fill in apparent accuracy
map.accuracy = accuracy(map);

