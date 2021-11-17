init;

% Set grid
seti.dim = 2;           % two dimensional problem
seti.rCD = 0.2;         % size of computational domain [-rCD,rCD)^dim
seti.nCD = 256;         % number of discretization points for each dimension of CD
seti = setGrid(seti);   % grid on CD is stored in seti.grid; grid on ROI is stored in seti.gridROI
                        % discretization points in each dimension are stored in seti.nROI
