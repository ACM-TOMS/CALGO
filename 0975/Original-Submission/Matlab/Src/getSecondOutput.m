% Get the second output from a function
%
% Example:
%
%  For function [x,y] = f(t)
%
%    val = getSecondOuput(@f,t)
%
%  returns y.
%
%  This is useful if one want an anonymous function that returns the second
%  ouptut value from f ie
%
%    g = @(t) getSecondOutput(@f,t).
%
% Stuart C. Hawkins - 6 November 2014

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


function val = getSecondOutput(fun,param)

[~,val] = fun(param);

