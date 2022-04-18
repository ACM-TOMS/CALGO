function C = gpopsuserfunFD(x);
%------------------------------------------------------------------%
% User function when numerical differentiation is used             %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao, Geoffrey T. Huntington, David   %
% Benson, Michael Patterson, Christopher Darby, & Camila Francolin %
%------------------------------------------------------------------%

global mysetup

%--------------------------------%
% Unscale the decision variables %
%--------------------------------%
y = x./mysetup.column_scales;
C = gpopsObjandCons(y);
%-----------------------------------------------------------------%
% Compute functions based on unscaled value of decision variables %
%-----------------------------------------------------------------%
C(2:mysetup.numnonlin+1) = mysetup.DF*C(2:mysetup.numnonlin+1);

