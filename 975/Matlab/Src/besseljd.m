% Derivative of Bessel function of the first kind.
%
%  z = besseljd(n,x) returns the derivative of the Bessel function of the
%  first kind of order n evaluated at points x. It is assumed that n is an
%  integer.
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


function val=besseljd(n,x)

% compute values using formula derived from (9.1.30) in Abramowitz and
% Stegun, Handboook of Mathematical functions.
val = besselj(n-1,x) - n.*besselj(n,x)./x;

% ensure the result is a column vector
val=val(:);