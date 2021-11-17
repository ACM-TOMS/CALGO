function cholmod_acm_toms_test
% cholmod_acm_toms_test exhaustive test of all CHOLMOD_ACM_TOMS packages
%
% Your current directory must be CHOLMOD_ACM_TOMS for this function to work.
% cholmod_acm_toms_install must be run prior to running this test.  Warning:
% this test takes a *** long **** time.
%
% Example:
%   cholmod_acm_toms_test
%
% See also cholmod_acm_toms_install, cholmod_acm_toms_demo.

% Copyright 2007, Tim Davis, University of Florida

npackages = 13 ;
h = waitbar (0, 'CHOLMOD ACM TOMS test:') ;
SuiteSparse = pwd ;

[v,pc] = getversion ;
if (v < 7)
    error ('cholmod_acm_toms_test requires MATLAB 7.0 or later') ;
end

% if at UF, ensure pre-installed UF Sparse Matrix Collection is used
uf = { '/cise/homes/davis/Install/UFget', 'd:/UFget', '/share/UFget' } ;
for k = 1:length(uf)
    if (exist (uf {k}, 'dir'))
        addpath (uf {k}) ;
        break ;
    end
end

try

    %---------------------------------------------------------------------------
    % COLAMD
    %---------------------------------------------------------------------------

    waitbar (3/(npackages+1), h, 'CHOLMOD ACM TOMS test: COLAMD') ;
    cd ([SuiteSparse '/COLAMD/MATLAB']) ;
    colamd_test ;

    %---------------------------------------------------------------------------
    % CCOLAMD
    %---------------------------------------------------------------------------

    waitbar (4/(npackages+1), h, 'CHOLMOD ACM TOMS test: CCOLAMD') ;
    cd ([SuiteSparse '/CCOLAMD/MATLAB']) ;
    ccolamd_test ;

    %---------------------------------------------------------------------------
    % CHOLMOD
    %---------------------------------------------------------------------------

    waitbar (6/(npackages+1), h, 'CHOLMOD ACM TOMS test: CHOLMOD') ;
    cd ([SuiteSparse '/CHOLMOD/MATLAB/Test']) ;
    cholmod_test ;

    %---------------------------------------------------------------------------
    % AMD, CAMD, UFget
    %---------------------------------------------------------------------------

    % no exhaustive tests; tested via other packages

catch

    %---------------------------------------------------------------------------
    % test failure
    %---------------------------------------------------------------------------

    cd (SuiteSparse) ;
    disp (lasterr) ;
    fprintf ('CHOLMOD ACM TOMS test: FAILED\n') ;
    return

end

%-------------------------------------------------------------------------------
% test OK
%-------------------------------------------------------------------------------

close (h) ;
fprintf ('CHOLMOD ACM TOMS test: OK\n') ;
cd (SuiteSparse) ;

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
