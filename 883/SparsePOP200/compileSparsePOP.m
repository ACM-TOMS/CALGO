% Compiling Mex Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a component of SparsePOP 
% Copyright (C) 2007 SparsePOP Project
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
eval('cd subPrograms/Mex');

MLVerStr = version;
MLVer = str2num(MLVerStr(1:3));

LIBfiles = ' conversion.cpp spvec.cpp polynomials.cpp sup.cpp ';
if ispc % Windows family create .obj files
        OBJfiles = ' conversion.obj spvec.obj polynomials.obj sup.obj ';
else
        OBJfiles = ' conversion.o spvec.o polynomials.o sup.o ';
end


if strcmp(computer, 'GLNXA64') && MLVer > 7.2
        MexFlags = ' -O -Dlinux=1 -largeArrayDims ';
elseif strcmp(computer, 'GLNX86') && MLVer > 7.2
        MexFlags = ' -O CXXFLAGS="-fPIC -ansi -D_GNU_SOURCE -pthread -m32 "';
elseif ispc
        MexFlags = ' -O -Dlinux=0 ';
else % Mac, Linux 32, or Solaris
        MexFlags = ' -O -Dlinux=0 ';
end

fprintf('Compiling Libraries...');
command = ['mex -c ' MexFlags LIBfiles];
eval(command);
fprintf('done\n');
fprintf('Generating mexconv1...');
command = ['mex ' MexFlags ' mexconv1.cpp'  OBJfiles ];
eval(command);
fprintf('done\n');
fprintf('Generating mexconv2...');
command = ['mex ' MexFlags ' mexconv2.cpp ' OBJfiles ];
eval(command);
fprintf('done\n');
eval('cd ../../');
fprintf('Compilation finished successfully.\n');

