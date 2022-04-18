
function Fadv = slt_interp(A,F,XD,YD)
% Semi-Lagrangian transport (interpolation) 
% of the field F at departure points (XD,YD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% A is the extended slt grid on which fields are defined 
%  in particular [A.X,A.Y] and [A.Xe A.Ye] are the meshgrid produced
%  arrival points and extended grid coordinate arrays
% F is nxm field of values to be advected from time level n
% XD,YD are the departure point values for which fields will be interpolated
% Output: 
% Fadv is the resulting field at the departure points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the values of F at the departure points
%[Xe,Ye] = meshgrid(x,y); % extended grid
Fe= slt_extend(A,F);
Fadv=interp2(A.Xe,A.Ye,Fe,XD,YD,'cubic'); % using high order interpolant
Fadv = Fadv';  %restore row dimension as lon
% end slt_interp
