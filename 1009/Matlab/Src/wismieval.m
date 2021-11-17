% Compute the recommended order for wavefunction expansion.
%
%   n = wismieval(s) returns the suggested expansion order for scatterers
%   with relative size s.
%
%   The relative size s = diameter / wavelength.
%
%   The suggested expansion order is based on the formula given on Page
%   1508 of W. J. Wiscombe, Applied Optics, Vol 19 No 9. 1 May 1980.
%
% Stuart C. Hawkins - 6 December 2018

% Copyright 2014-2019 Stuart C. Hawkins
% 	
% This file is part of MIESOLVER.
% 
% MIESOLVER is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MIESOLVER is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MIESOLVER.  If not, see <http://www.gnu.org/licenses/>.


function nmie=wismieval(relsize)

%---------------------------
% set wavenumber etc
%---------------------------

% This is the strategy from
% Improved Mie Scattering Algorithms,
% W. J. Wiscombe
% Applied Optics, Vol 19 No 9. 1 May 1980.
% Page 1508.

% compute size parameter
x=pi*relsize;

% set default
if x<=8
    nmie=floor(x)+1+ceil(4*x^(1/3));
elseif x<4200
    nmie=floor(x)+2+ceil(4.05*x^(1/3));
else
    nmie=floor(x)+2+ceil(4*x^(1/3));
end
