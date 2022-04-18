% Enumerated type for holding the kind of boundary condition.
%
% See http://stackoverflow.com/questions/1389042/how-do-i-create-enumerated-types-in-matlab

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


classdef (Sealed) boundary_condition
    properties (Constant)
        soft = 0;
        hard = 1;
        robin = 2;
        transmission = 3;
    end

    methods (Access = private)    % private so that you can't instantiate
        function out = boundary_condition
        end
    end
end