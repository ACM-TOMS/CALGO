% Suggested order for wavefunction expansion
%
%  n = suggestedorder(k,r) returns the suggested order n for a regular 
%  wavefunction expansion of a plane wave with wavenumber k valid in a ball
%  of radius r.
%
% The suggested order is computed from a formula on p1508 of Improved Mie
% Scattering Algorithms, W. J. Wiscombe, Applied Optics Vol 19 No 9, 1980.
%
% See also: regularwavefunctionexpansion.
%
% Stuart C. Hawkins - 12 November 2014

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


function val = suggestedorder(kwave,radius)

% This is the strategy from
% Improved Mie Scattering Algorithms,
% W. J. Wiscombe
% Applied Optics, Vol 19 No 9. 1 May 1980.
% Page 1508.

% compute size parameter
x=kwave*radius;

% set default
if x<=8
    val=floor(x)+1+ceil(4*x^(1/3));
elseif x<4200
    val=floor(x)+2+ceil(4.05*x^(1/3));
else
    val=floor(x)+2+ceil(4*x^(1/3));
end
