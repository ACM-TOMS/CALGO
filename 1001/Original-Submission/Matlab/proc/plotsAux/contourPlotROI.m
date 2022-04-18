%% contourPlotROI
% Restricts the contrast to ROI and plots 3D array q on ROI as contour.
%
%% Syntax
%
%   qROI = contourPlotROI(q, seti, part, levelPlotValue)
%
%% Description
%
% |qROI = contourPlotROI(q, seti, part, levelPlotValue)| restricts the
% contrast |q| to ROI and plots isosurface data with isosurface value
% |levelPlotValue| of the real or imaginary part of the contrast 
% (|part| is 'real' or 'imag').
%
%% Input Arguments
%
% * q    : data (vector or array
% * seti : structural array
% * part : real or imaginary part: 'real' or 'imag'
% * levelPlotValue : level in [0,1] to plot levelset (times max(abs(q))
%                    (This is an optional argument and set to 0.3 in case of missing).
%
% The structural array |seti| has to contain the fields:
%
% * seti.ROImask, 
% * seti.ballMaskROI, 
% * seti.gridROI.
%
% For details of this fields see <setGrid.html>.
%
%% Output Arguments
%
% * qROI: restriction of q to ROI
%
% Plot the isosurface data with isosurface value |levelPlotValue|.
%
%% See Also:
%
% * <plotAndSaveFigures.html>
% * <setGrid.html>
%
% * MathWorks. MATLAB documentation.
%   <https://de.mathworks.com/help/matlab/ref/contour.html>.
%   Accessed: 2016-11-03.
% *  MathWorks. MATLAB documentation.
%   <https://de.mathworks.com/help/matlab/ref/isosurface.html>.
%   Accessed: 2016-11-03.
%
%% Code
function qROI = contourPlotROI(q, seti, part, levelPlotValue)

% part: 'real' or 'imag'
%if ~exist(part,'var')
%    part = 'real';
%end

%%
% *Resctrict CD to ROI*

if (length(size(q))==2)&&(numel(q)==max(size(q))) % data is vector 
    if numel(q)==numel(seti.ROImask)
        q = reshape(q, size(seti.ROImask));
        qROI = restrictCDtoROI(q, seti.ROImask);
    elseif numel(q)==numel(seti.ballMaskROI)
        qROI = reshape(q, size(seti.ballMaskROI));
    else
        fprintf('Error in contourPlotROI.m: Do not find suitable size of vector to transform vector into array\n');
        qROI=q;
    end
elseif length(size(q))==2 % data is matrix 
    if numel(q)==numel(seti.ROImask)
        qROI = restrictCDtoROI(q, seti.ROImask);
    elseif numel(q)==numel(seti.ballMaskROI)
        qROI = reshape(q, size(seti.ballMaskROI));
    else
        fprintf('Error in contourPlotROI.m: Do not find suitable size of 2D array\n');
        qROI=q;
    end
elseif length(size(q))==3
    if numel(q)==numel(seti.ROImask)
        qROI = restrictCDtoROI(q, seti.ROImask);
    elseif numel(q)==numel(seti.ballMaskROI)
        qROI=q;
    else
        fprintf('Error in contourPlotROI.m: Do not find suitable size of 3D array\n');
        qROI=q;
    end
end
% qROI should now be the restriction of q to ROI and transformed into seti.dim-array

%%
% *Set levelPlotValue*

if ~exist('levelPlotValue', 'var')
    levelPlotValue = 0.3; % default
end

%%
% *Real or imaginary part*

% plot real part
% qROI=real(qROI); 

switch part
    case 'real'
        qROI = real(qROI);
    case 'imag'
        qROI = imag(qROI);
    otherwise
        error('You have to choose real or imag part.')
end

%%
% *Isosurface*

p = patch(isosurface(reshape(seti.gridROI(1,:),size(seti.ballMaskROI)), ...
    reshape(seti.gridROI(2,:),size(seti.ballMaskROI)), ...
    reshape(seti.gridROI(3,:),size(seti.ballMaskROI)), qROI , levelPlotValue));
isonormals(reshape(seti.gridROI(1,:),size(seti.ballMaskROI)), ...
    reshape(seti.gridROI(2,:),size(seti.ballMaskROI)), ...
    reshape(seti.gridROI(3,:),size(seti.ballMaskROI)), qROI, p);
set(p,'FaceColor','red','EdgeColor','none');
daspect([1,1,1])
view(3); 
axis([min(seti.gridROI(1,:)) max(seti.gridROI(1,:)) min(seti.gridROI(2,:)) ...
      max(seti.gridROI(2,:)) min(seti.gridROI(3,:)) max(seti.gridROI(3,:))]);
grid on
camlight 
daspect([1 1 1])
camlight(-45,50)
hold on
xlabel('x');
ylabel('y');
zlabel('z');
%title(sprintf('reconstructed contrast (%s part)',part));

X_P = get(p,'XData');
Y_P = get(p,'YData');
Z_P = get(p,'ZData');

s1 = patch(max(seti.gridROI(1,:))*ones(size(Y_P)),Y_P,Z_P,[0.75,0.75,0.75]);
s2 = patch(X_P,max(seti.gridROI(2,:))*ones(size(X_P)),Z_P,[0.75,0.75,0.75]);
s3 = patch(X_P,Y_P,min(seti.gridROI(3,:))*ones(size(X_P)),[0.25,0.25,0.25]);

set(s1,'EdgeColor','none');
set(s2,'EdgeColor','none');
set(s3,'EdgeColor','none');

hold off

end
