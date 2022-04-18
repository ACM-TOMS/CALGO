function SuiteSparse_install (do_demo)
%SuiteSparse_install: compiles and installs all of SuiteSparse
% (SuiteSparseQR subset: SuiteSparseQR and all its dependent packages)
% A Suite of Sparse matrix packages, authored or co-authored by Tim Davis, Univ.
% Florida. You must be in the same directory as SuiteSparse_install to use this.
%
% Packages in SuiteSparse (SuiteSparseQR subset):
%
% CHOLMOD        sparse Cholesky factorization, and many other operations
% AMD            sparse symmetric approximate minimum degree ordering
% COLAMD         sparse column approximate minimum degree ordering
% CAMD           constrained AMD
% CCOLAMD        constrained COLAMD
% SuiteSparseQR  sparse QR factorization
%
% Except where noted, all packages work on MATLAB 6.1 or later.  They have not
% been tested on earlier versions, but they might work there.  Please let me
% know if you try SuiteSparse on MATLAB 6.0 or earlier, whether it works or not.
%
% Example:
%    SuiteSparse_install
%    help SuiteSparse      % for more details
%
% See also AMD, COLAMD, CAMD, CCOLAMD, CHOLMOD, SuiteSparse, SPQR,
%   PATHTOOL, PATH.

% Copyright 1990-2008, Timothy A. Davis.
% http://www.cise.ufl.edu/research/sparse
% In collaboration with Patrick Amestoy, Yanqing Chen, Iain Duff, John Gilbert,
% Steve Hadfield, Bill Hager, Stefan Larimore, Esmond Ng, Eka Palamadai, and
% Siva Rajamanickam.

paths = { } ;
SuiteSparse = pwd ;

% determine the MATLAB version (6.1, 6.5, 7.0, ...)
v = getversion ;
pc = ispc ;

% check if METIS 4.0.1 is present where it's supposed to be
have_metis = exist ('metis-4.0', 'dir') ;
if (~have_metis)
    fprintf ('Some codes optionally use METIS 4.0.1.  Download\n') ;
    fprintf ('it from http://glaros.dtc.umn.edu/gkhome/views/metis\n');
    fprintf ('and place the metis-4.0 directory in this directory.\n') ;
    input ('or hit enter to continue without METIS: ', 's') ;
    fprintf ('Now compiling without METIS...\n\n') ;
end

% print the introduction
help SuiteSparse_install

fprintf ('MATLAB version %g (%s)\n', v, version) ;

% add SuiteSparse to the path
fprintf ('\nPlease wait while SuiteSparse is compiled and installed...\n') ;
paths = add_to_path (paths, SuiteSparse) ;

% compile and install CHOLMOD
try
    % determine whether or not to compile CHOLMOD
    cd ([SuiteSparse '/CHOLMOD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    if (have_metis)
       cholmod_make
    else
       cholmod_make ('no metis') ;
    end
catch                                                                       %#ok
    fprintf ('CHOLMOD not installed\n') ;
end

% compile and install AMD
try
    cd ([SuiteSparse '/AMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    amd_make
catch                                                                       %#ok
    fprintf ('AMD not installed\n') ;
end

% compile and install COLAMD
try
    cd ([SuiteSparse '/COLAMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    colamd_make
catch                                                                       %#ok
    fprintf ('COLAMD not installed\n') ;
end

% compile and install CCOLAMD
try
    cd ([SuiteSparse '/CCOLAMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    ccolamd_make
catch                                                                       %#ok
    fprintf ('CCOLAMD not installed\n') ;
end

% compile and install CAMD
try
    cd ([SuiteSparse '/CAMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    camd_make
catch                                                                       %#ok
    fprintf ('CAMD not installed\n') ;
end

% compile and install SuiteSparseQR
try
    if (pc)
        fprintf ('Note that SuiteSparseQR will not compile with the lcc\n') ;
        fprintf ('compiler provided with MATLAB on Windows\n') ;
    end
    cd ([SuiteSparse '/SPQR/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    if (have_metis)
       spqr_make
    else
       spqr_make ('no metis') ;
    end
catch                                                                       %#ok
    disp (lasterr) ;                                                        %#ok
    fprintf ('SuiteSparseQR not installed\n') ;
end

% post-install wrapup

cd (SuiteSparse)
fprintf ('SuiteSparse is now installed.\n') ;

if (nargin < 1)
    % ask if demo should be run
    y = input ('Hit enter to run the SuiteSparse demo (or "n" to quit): ', 's') ;
    if (isempty (y))
        y = 'y' ;
    end
    do_demo = (y (1) ~= 'n') ;
end
if (do_demo)
    try
	SuiteSparse_demo ;
    catch                                                                   %#ok
	fprintf ('SuiteSparse demo failed\n') ;
    end
end

fprintf ('\nSuiteSparse installation is complete.  The following paths\n') ;
fprintf ('have been added for this session.  Use pathtool to add them\n') ;
fprintf ('permanently.  If you cannot save the new path because of file\n');
fprintf ('permissions, then add these commands to your startup.m file.\n') ;
fprintf ('Type "doc startup" and "doc pathtool" for more information.\n\n') ;
for k = 1:length (paths)
    fprintf ('addpath %s\n', paths {k}) ;
end
cd (SuiteSparse)

fprintf ('\nSuiteSparse for MATLAB %g installation complete\n', getversion) ;

%-------------------------------------------------------------------------------
function paths = add_to_path (paths, newpath)
% add a path
addpath (newpath) ;
paths = [paths { newpath } ] ;						    %#ok

%-------------------------------------------------------------------------------
function v = getversion
% determine the MATLAB version, and return it as a double.
v = sscanf (version, '%d.%d.%d') ;
v = 10.^(0:-1:-(length(v)-1)) * v ;
