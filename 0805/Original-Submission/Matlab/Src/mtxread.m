function [A, m, n] = mtxread(infile);

%MTXREAD Read matrix A from a file in matrix market format.
%
%   A = mtxread(FILENAME) reads a MatrixMarket formatted matrix from the
%   file FILENAME. 
%
%SDDPACK: Software for the Semidiscrete Decomposition.
%Copyright (c) 1999 Tamara G. Kolda and Dianne P. O'Leary. 

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free 
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.  
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.  
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc., 59
% Temple Place - Suite 330, Boston, MA 02111-1307, USA.  

fid = fopen(infile, 'rt');

if (fid == -1)
  error('Error opening file.');
end

line = fgets(fid);
while line(1) == '%'
  line = fgets(fid);
end

[data, cnt] = sscanf(line,'%d');
m = data(1);
n = data(2);
nnzs = data(3);

[data, cnt] = fscanf(fid, '%d %d %e', [3, inf]);

I = data(1,:);
J = data(2,:);
S = data(3,:);
A = sparse(I, J, S, m, n);

fclose(fid);

return;
