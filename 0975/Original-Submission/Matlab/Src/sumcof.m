% Evaluate wavefunction expansion series.
%
%   z = sumcof(x,x0,k,c,'H') returns the values z of the radiating 
%   wavefunction expansion with coefficients c, centre x0 and wavenumber k 
%   at points x.
%
%   z = sumcof(x,x0,k,c,'J') returns the values z of the regular 
%   wavefunction expansion with coefficients c, centre x0 and wavenumber k 
%   at points x.
%
%   z = sumcof(x,x0,k,c,'F') returns the values z of the far field of the 
%   radiating wavefunction expansion with coefficients c, centre x0 and 
%   wavenumber k at points abs(x) on the unit circle.
%
% Note: in the above vectors in the plane are represented by
% complex numbers.
%
% Stuart C. Hawkins - 13 January 2015

% Copyright 2014, 2015 Stuart C. Hawkins and M. Ganesh.
% 	
% This file is part of TMATROM.
% 
% TMATROM is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TMATROM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TMATROM.  If not, see <http://www.gnu.org/licenses/>.


function val = sumcof(points,centre,kwave,cof,type)

%-------------------------------------------------
% setup
%-------------------------------------------------

% make sure the coefficient vector is a column vector
cof=cof(:);

% determine the maximum order from the length of the coefficient vector
nmax=0.5*(length(cof)-1);

% create a vector of indexes.... helps to vectorize the computation
n=-nmax:nmax;

%-------------------------------------------------
% turn points into a vector
%-------------------------------------------------

% get the shape of points so we can restore it later
[np,mp]=size(points);

% reshape into a vector
if strcmp(type,'F')
    p=reshape(points,np*mp,1);
else
    p=reshape(points-centre,np*mp,1);
end

%-------------------------------------------------
% compute the field
%-------------------------------------------------

% convert to polar coordinates
theta=angle(p);
rad=abs(p);

% make a matrix from n and rad
[nd,rd]=meshgrid(n,kwave*rad);

% get Bessel/Hankel/far field values as appropriate
if strcmp(type,'J')
    
    bess=besselj(abs(nd),rd);
    
elseif strcmp(type,'H')
    
    bess=besselh(abs(nd),rd);
    
elseif strcmp(type,'F')
    
    bess=sqrt(1/(pi*kwave))*(-1i).^abs(nd)*(1-1i);   
    
end

% compute the angular part
Y=exp(1i*theta*n);

% adjust the farfield if the centre of the scatterer is not
% the origin
if strcmp(type,'F') && centre~=0
    phase = exp(-1i*kwave*real(exp(-1i*theta)*centre));
    Y = diag(phase)*Y;
end

% put it together
M = bess.*Y;

%-------------------------------------------------
% make the return value the same shape as the original
% array of points
%-------------------------------------------------

% compute the sum of the wavefunctions and reshape
val=reshape(M*cof,np,mp);