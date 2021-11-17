%% inseti: example
% Example of input parameters to run |start| not with default parameters.

%% Description
% If you want to run |start| not with the default parameters you can set
% the fields in struct |seti| in this or another file in the folder
% |inseti| and refer to this file by Parameter |inseti|.
%
% *Attention in usage of the framework*
%
% * |start| calls |setInput|, in which *pre-defined variables* (except closed test and inseti) *are deleted*.
% * Therefore you *must define* the struct *|seti| in* such a *file* (and not in terminal).
% * The name of a *new file in |inseti|* must differ from existing functions
% (because the file in |inseti| is called by |eval|).
%
%% Syntax and Example
%   inseti = 'example';
%   start
%
%% Input Parameters of struct |seti|
% Note that this are *some*, but *not all* input parameters.

% -- dirname: suffix
seti.dirSuffix = '_test';   % Suffix of the folder created in folder |output| to store the results

% -- setGrid
seti.dim = 2;       % dimension of scattering problem (2 or 3)
                    % Make sure to choose a suitable contrast in seti.contrast.
seti.rCD = 0.2;     % size of computational domain |[-rCD,rCD]^dim|
seti.nCD = 256;     % number of discretization points in each dimension

% -- setKernel
seti.k = 250;               % wave number
seti.model = 'helmholtz';   % input helmholtz will choose automatically 
                            % helmholtz2D or helmholtz3D dependent on the dim

% -- setContrast
seti.contrast = 'cornerBallSparse2D';
% Set the name of contrast function as string (files in folder |incontrasts|).
% Make sure to choose the correct dimension of the problem in |seti.dim|.
% A list of available contrasts is in setContrast.

% -- expSetup (set experimental set-up)
seti.incPntsType = 'circle';    % Geometry: type of transmitters
seti.measPntsType = 'circle';   % Geometry: type of receivers
seti.incNb = 35;                % Geometry: number of transmitters
seti.measNb = 35;               % Geometry: number of receivers
seti.radSrc = 5;                % Geometry: radius of sphere containing transmitters
seti.radMeas = 5;               % Geometry: radius of sphere containing receivers
% Fore more details of geometry type circle and further geometries like square and borehole, see |pntsGeometry|.
seti.incType = 'pointSource';   % Type of incident fields: 'pointSource' or 'planeWave'
seti.measType = 'nearField';    % Type of measurements: 'nearField' or 'farField'

% -- reconstruction
seti.invNo = 6;                 % Do not change in published code (no other option).
seti.delta = 0.01;              % Relative noise level (for noisy simulated data and assumed in reconstruction process to stop by discrepancy principle)
seti.physBounds = [-1,3,0,3];   % Bounds for real/imaginary part of contrast: [reMin, reMax, imMin, imMax]
seti.alpha = 500;               % Regularization parameter of sparse penalty term
seti.beta = 1E-5;               % Regularization parameter of total variation penalty term
seti.useDis = 1;                % 1, then discrepancy principle is used to stop reconstruction process (outer iteration)
seti.tau = 2.5;                 % Discrepancy principle stops at parameter \tau \delta
seti.nOut = 30;                 % Maximal number of reconstruction steps (outer iteration)

% -- reconstruction with PDA
seti.pdaN = 50;                 % number of inner iteration steps (PDA)
seti.pdaStepsize = 'fix';       % how primal and dual stepsizes in pda are choosen
                                % 'fix' (tau = sigma = 1/L).
                                % In our experience you should not change it.

% -- figures and files
seti.plotFreq = 1;              % Frequency to plot outer iteration (0: no plot)
seti.plotPublish = 1;           % 1: plot design for publication (e.g. without title)
                                % see plot2DstylePublish.m and plot3DstylePublish.m
seti.usecbarlim = 1;            % 1: set limits of colorbar manually (0: automatically)
seti.cbarlim = [-0.2, 1.4];     % limits vor colorbars [cbarmin, cbarmax]

seti.savepng = 1;               % 0 or 1, save figures as *.png (default: 1)
seti.saveepsc = 0;              % 0 or 1, save figures as colored *.eps (default: 0)
seti.savefig = 0;               % 0 or 1, save figures as *.fig (default: 0)

seti.savedata = 1;              % 0 or 1, saves relative discrepancy, error 
                                % and difference of computed contrast to its predecessor
                                % as save_dis.mat, save_err.mat and
                                % save_dif.mat in folder output.

%% More About
%
% * |rCD|: In comparison to a similar Parameter $R$ in the Manuscript, see 
% [1], we have the relation: $R$ = |rCD/2| = 0.1.

%% References
% [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.

%% See Also
%
% * <start.html>
% * <setContrast.html>
% * <pntsGeometry.html>
