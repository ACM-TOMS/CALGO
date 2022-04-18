function SuiteSparse_install (do_demo)
%SuiteSparse_install: compiles and installs all of SuiteSparse
% A Suite of Sparse matrix packages, authored or co-authored by Tim Davis.
%
% Packages in SuiteSparse (sparse QR subset):
%
% CHOLMOD        sparse Cholesky factorization, and many other operations
% AMD            sparse symmetric approximate minimum degree ordering
% COLAMD         sparse column approximate minimum degree ordering
% CAMD           constrained AMD
% CCOLAMD        constrained COLAMD
% UFget          interface to SuiteSparse Matrix Collection
% SuiteSparseQR  sparse QR factorization
%
% Example:
%    SuiteSparse_install        % compile and prompt to run each package's demo
%    SuiteSparse_install(0)     % compile but do not run the demo
%    SuiteSparse_install(1)     % compile and run the demos with no prompts
%    help SuiteSparse           % for more details
%
% See also AMD, COLAMD, CAMD, CCOLAMD, CHOLMOD, UMFPACK, CSPARSE, CXSPARSE,
%      UFget, RBio, UFcollection, KLU, BTF, MESHND, SSMULT, LINFACTOR, SPOK,
%      SPQR_RANK, SuiteSparse, SPQR, PATHTOOL, PATH, FACTORIZE, SPARSEINV.
%
% Copyright 1990-2016, Timothy A. Davis, http://www.suitesparse.com.
% In collaboration with Patrick Amestoy, Yanqing Chen, Iain Duff, John Gilbert,
% Steve Hadfield, William Hager, Stefan Larimore, Leslie Foster, Eka Palamadai
% Natarajan, Esmond Ng, Siva Rajamanickam, Nuri Yeralan, Sanjay Ranka,
% and Wissam Sid-Lakhdar.

%-------------------------------------------------------------------------------
% initializations
%-------------------------------------------------------------------------------

paths = { } ;
SuiteSparse = pwd ;

% determine the MATLAB version (6.1, 6.5, 7.0, ...)
v = version ;
pc = ispc ;

% print the introduction
help SuiteSparse_install

fprintf ('\nInstalling SuiteSparse for MATLAB version %s\n\n', v) ;

% add SuiteSparse to the path
paths = add_to_path (paths, SuiteSparse) ;

%-------------------------------------------------------------------------------
% compile and install the packages
%-------------------------------------------------------------------------------

% compile and install CHOLMOD
try
    paths = add_to_path (paths, [SuiteSparse '/CHOLMOD/MATLAB']) ;
    cholmod_make ;
catch me
    disp (me.message) ;
    fprintf ('CHOLMOD not installed\n') ;
end

% compile and install AMD
try
    paths = add_to_path (paths, [SuiteSparse '/AMD/MATLAB']) ;
    amd_make ;
catch me
    disp (me.message) ;
    fprintf ('AMD not installed\n') ;
end

% compile and install COLAMD
try
    paths = add_to_path (paths, [SuiteSparse '/COLAMD/MATLAB']) ;
    colamd_make ;
catch me
    disp (me.message) ;
    fprintf ('COLAMD not installed\n') ;
end

% compile and install CCOLAMD
try
    paths = add_to_path (paths, [SuiteSparse '/CCOLAMD/MATLAB']) ;
    ccolamd_make ;
catch me
    disp (me.message) ;
    fprintf ('CCOLAMD not installed\n') ;
end

% compile and install CAMD
try
    paths = add_to_path (paths, [SuiteSparse '/CAMD/MATLAB']) ;
    camd_make ;
catch me
    disp (me.message) ;
    fprintf ('CAMD not installed\n') ;
end

% install UFget, unless it's already in the path
try
    % if this fails, then UFget is not yet installed
    index = UFget ;
    fprintf ('UFget already installed:\n') ;
    which UFget
catch
    index = [ ] ;
end
if (isempty (index))
    % UFget is not installed.  Use SuiteSparse/UFget
    fprintf ('Installing SuiteSparse/UFget\n') ;
    try
        paths = add_to_path (paths, [SuiteSparse '/UFget']) ;
    catch me
        disp (me.message) ;
        fprintf ('UFget not installed\n') ;
    end
end

% compile and install SuiteSparseQR
try
    if (pc)
        fprintf ('Note that SuiteSparseQR will not compile with the lcc\n') ;
        fprintf ('compiler provided with MATLAB on Windows\n') ;
    end
    paths = add_to_path (paths, [SuiteSparse '/SPQR/MATLAB']) ;
    spqr_make ;
catch me
    disp (me.message) ;
    fprintf ('SuiteSparseQR not installed\n') ;
end

%-------------------------------------------------------------------------------
% post-install wrapup
%-------------------------------------------------------------------------------

cd (SuiteSparse)
fprintf ('SuiteSparse is now installed.\n\n') ;

% run the demo, if requested
if (nargin < 1)
    % ask if demo should be run
    y = input ('Hit enter to run the SuiteSparse demo (or "n" to quit): ', 's');
    if (isempty (y))
        y = 'y' ;
    end
    do_demo = (y (1) ~= 'n') ;
    do_pause = true ;
else
    % run the demo without pausing
    do_pause = false ;
end
if (do_demo)
    try
	SuiteSparse_demo ([ ], do_pause) ;
    catch me
        disp (me.message) ;
	fprintf ('SuiteSparse demo failed\n') ;
    end
end

% print the list of new directories added to the path
fprintf ('\nSuiteSparse installation is complete.  The following paths\n') ;
fprintf ('have been added for this session.  Use pathtool to add them\n') ;
fprintf ('permanently.  If you cannot save the new path because of file\n');
fprintf ('permissions, then add these commands to your startup.m file.\n') ;
fprintf ('Type "doc startup" and "doc pathtool" for more information.\n\n') ;
for k = 1:length (paths)
    fprintf ('addpath %s\n', paths {k}) ;
end
cd (SuiteSparse)

fprintf ('\nSuiteSparse for MATLAB %s installation complete\n', v) ;

%-------------------------------------------------------------------------------
function paths = add_to_path (paths, newpath)
% add a path
cd (newpath) ;
addpath (newpath) ;
paths = [paths { newpath } ] ;
