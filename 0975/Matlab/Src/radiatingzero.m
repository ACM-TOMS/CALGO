% Radiating wavefunction with zero coefficients.
%
%  u = radiatingzero(n,x,k) returns a radiating wavefunction expansion 
%  object u with order n, expansion origin x, wavenumber k and zero 
%  coefficient vector.
%
% Note: in the above vectors in the plane are represented by
% complex numbers.
%
% See also: regularzero, radiatingwavefunctionexpansion.
%
% Stuart C. Hawkins - 4 March 2015

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


function val = radiatingzero(n,x0,k)

val = radiatingwavefunctionexpansion(n,x0,k,zeros(2*n+1,1));