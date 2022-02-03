
function Fe = slt_extend(A,F)
% Extend the field F for SLT interpolation on the Xe,Ye grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% A the slt_grid extended grid
% A.X,Y are the grid values on which fields are defined produced by meshgrid
% F is nxm array of interior values to be extended
% Output: 
% Fe is the extend field over the extended grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = A.ni;  % interior index ranges
ny = A.nj;
nxpt = A.nxpt;
Fe = zeros(A.nie,A.nje); 
Fe( (nxpt+1):(nx+nxpt), (nxpt+1):(ny+nxpt) ) = F;  % fill the interior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -------------------------------------
%  |                                   |
%  |       C                D          |
%  |  +++++++++++++++N++++++++++++++   +1.0
%  |  +                            +   |
%  |  +                            +   |
%  |  +            interior        + periodic
%  |  +                            +   |
%  |  +                            +   |
%  |  +++++++++++++++S++++++++++++++   -1.0 (= sin lat)
%  | 0.0   A                B      2*pi
%  |                                   |
%  |                                   |
%  -------------------------------------
%     +nxpt+1                      +nxpt+nx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill periodic halos (West-East or Left-Right)
for i = 1:nxpt
	Fe(i,:) = Fe(nx+i,:);
	Fe(nxpt+nx+i,:) = Fe(nxpt+i,:);
end
% Fill pole extensions (South-North or Bottom-Top) 
%  fill regions are A(South-West), B(South-East), C(North-West), D(North-East)
%  (Better check this carefully.)
nmid = A.nie/2;
nend = A.nie;
for j = 1:nxpt
	Fe(1:nmid,j)         = Fe((nmid+1):nend,2*nxpt-j+1);   % region A
	Fe((nmid+1):nend,j)  = Fe(1:nmid,2*nxpt-j+1);          % region B
	Fe(1:nmid,ny+nxpt+j) = Fe((nmid+1):nend,ny+nxpt-j+1);  % region C
	Fe((nmid+1):nend,ny+nxpt+j) = Fe(1:nmid,ny+nxpt-j+1);  % region D
end
% reshape for interpolation routines
Fe = Fe';
%  end function slt_extend
