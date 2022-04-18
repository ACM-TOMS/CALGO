
function val = get(A,prop)
%GET  gets a value of the slt_grid A
%  Usage:  val = get(A,prop)
%    Once the extended grid is established, with index mappings for
%    halo updates, then particle tracking and interpolation at departure
%    points can be done from input fields.
% Input:  A  slt_grid
%         prop  - property, i.e.  
% Output: val
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The extended grid data structure contains
%    ni, nj           the number of internal lons and lats
%    nie, nje         the number of total lons and lats in extend slt_grid
%    xge,yge  -   the extended grid in the x and y directions
%    nxpt 	-the number of extended points in the halo
%    X,Y   	-meshgrid output of the interior grid
%    Xe,Ye   	-meshgrid output of the extended grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isa(A,'slt_grid'))
	switch prop
	case 'xge'
	  val = A.xge;
	case 'yge'
	  val = A.yge;
	case 'ni'
	  val = A.ni;
	case 'nj'
	  val = A.nj;
	case 'nie'
	  val = A.nie;
	case 'nxpt'
	  val = A.nxpt;
	case 'nypt'
	  val = A.nypt;
	case 'X'
	  val = A.X;
	case 'Y'
	  val = A.Y;
	case 'Xe'
	  val = A.Xe;
	case 'Ye'
	  val = A.Ye;
        otherwise
           error('Get slt_grid method called with unknown property')
 	end
else
  error('Call to get slt_grid class method with wrong type')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end get
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

