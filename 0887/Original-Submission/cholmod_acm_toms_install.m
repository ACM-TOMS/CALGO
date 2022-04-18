function cholmod_acm_toms_install (do_demo)
%cholmod_acm_toms_install: compiles and installs CHOLMOD and its orderings
% You must be in the same directory as cholmod_acm_toms_install to use this.
%
% Packages in CHOLMOD_ACM_TOMS:
%
% CHOLMOD        sparse Cholesky factorization, and many other operations
% AMD            sparse symmetric approximate minimum degree ordering
% COLAMD         sparse column approximate minimum degree ordering
% CAMD           constrained AMD
% CCOLAMD        constrained COLAMD
% UFget          interface to UF Sparse Matrix Collection (MATLAB 7.0 or later)
%
% Except where noted, all packages work on MATLAB 6.1 or later.  They have not
% been tested on earlier versions, but they might work there.  Please let me
% know if you try CHOLMOD_ACM_TOMS on MATLAB 6.0 or earlier, whether it works
% or not.
%
% Example:
%    cholmod_acm_toms_install
%
% See also AMD, COLAMD, CAMD, CCOLAMD, CHOLMOD, UFget, PATHTOOL, PATH.

% Copyright 1990-2007, Timothy A. Davis.
% http://www.cise.ufl.edu/research/sparse
% In collaboration with Patrick Amestoy, Yanqing Chen, Iain Duff, John Gilbert,
% Steve Hadfield, Bill Hager, Stefan Larimore, Esmond Ng, Eka Palamadai, and
% Siva Rajamanickam.

paths = { } ;

% determine the MATLAB version (6.1, 6.5, 7.0, ...)
[v,pc] = getversion ;

% check if METIS 4.0.1 is present where it's supposed to be
have_metis = exist ('metis-4.0', 'dir') ;
if (~have_metis)
    fprintf ('CHOLMOD optionally uses METIS 4.0.1.  Download it\n') ;
    fprintf ('from http://glaros.dtc.umn.edu/gkhome/views/metis\n');
    fprintf ('and place the metis-4.0 directory in this directory.\n') ;
    fprintf ('Now compiling without METIS...\n\n') ;
end

% print the introduction
help cholmod_acm_toms_install

fprintf ('MATLAB version %.1f (%s)\n', v, version) ;

% add SuiteSparse to the path
SuiteSparse = pwd ;
fprintf ('\nPlease wait while CHOLMOD_ACM_TOMS is compiled and installed...\n');
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
catch
    fprintf ('CHOLMOD not installed\n') ;
end

% compile and install AMD
try
    cd ([SuiteSparse '/AMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    amd_make
catch
    fprintf ('AMD not installed\n') ;
end

% compile and install COLAMD
try
    cd ([SuiteSparse '/COLAMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    colamd_make
catch
    fprintf ('COLAMD not installed\n') ;
end

% compile and install CCOLAMD
try
    cd ([SuiteSparse '/CCOLAMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    ccolamd_make
catch
    fprintf ('CCOLAMD not installed\n') ;
end

% compile and install CAMD
try
    cd ([SuiteSparse '/CAMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    camd_make
catch
    fprintf ('CAMD not installed\n') ;
end

% compile and install CXSparse and UFget
try
    if (v >= 7.0)
	paths = add_to_path (paths, [SuiteSparse '/UFget']) ;
	fprintf ('UFget installed successfully\n') ;
    else
	fprintf ('UFget skipped; requires MATLAB 7.0 or later\n') ;
    end
catch
    fprintf ('UFget not installed\n') ;
end

% post-install wrapup
cd (SuiteSparse)
fprintf ('CHOLMOD_ACM_TOMS is now installed.\n') ;

if (nargin < 1)
    % ask if demo should be run
    y = input ('Hit enter to run the demo (or "n" to quit): ', 's') ;
    if (isempty (y))
        y = 'y' ;
    end
    do_demo = (y (1) ~= 'n') ;
end
if (do_demo)
    try
	cholmod_acm_toms_demo ;
    catch
	fprintf ('CHOLMOD demo failed\n') ;
    end
end

fprintf ('\nCHOLMOD_ACM_TOMS installation is complete.  The following paths\n');
fprintf ('have been added for this session.  Use pathtool to add them\n') ;
fprintf ('permanently.  If you cannot save the new path because of file\n');
fprintf ('permissions, then add these commands to your startup.m file.\n') ;
fprintf ('Type "doc startup" and "doc pathtool" for more information.\n\n') ;
for k = 1:length (paths)
    fprintf ('addpath %s\n', paths {k}) ;
end
cd (SuiteSparse)


%-------------------------------------------------------------------------------
function paths = add_to_path (paths, newpath)
% add a path
addpath (newpath) ;
paths = [paths { newpath } ] ;						    %#ok

%-------------------------------------------------------------------------------
function [v,pc] = getversion
% determine the MATLAB version, and return it as a double.
% only the primary and secondary version numbers are kept.
% MATLAB 7.0.4 becomes 7.0, version 6.5.2 becomes 6.5, etc.
v = version ;
t = find (v == '.') ;
if (length (t) > 1)
    v = v (1:(t(2)-1)) ;
end
v = str2double (v) ;
try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
catch
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
end
