function cholmod_acm_toms_demo
%CHOLMOD_ACM_TOMS_DEMO a demo of all packages in CHOLMOD_ACM_TOMS
%
% Example:
%   cholmod_acm_toms_demo
%
% See also cholmod, amd, camd, colamd, ccolamd

% Copyright (c) Timothy A. Davis, Univ. of Florida

input ('Hit enter to run the CHOLMOD demo: ') ;
try
    cholmod_demo
catch
    disp (lasterr) ;
    fprintf ('CHOLMOD demo failed\n' )
end

input ('Hit enter to run the CHOLMOD graph partitioning demo: ') ;
try
    graph_demo
catch
    disp (lasterr) ;
    fprintf ('graph_demo failed, probably because METIS not installed\n') ;
end

input ('Hit enter to run the AMD demo: ') ;
try
    amd_demo
catch
    disp (lasterr) ;
    fprintf ('AMD demo failed\n' )
end

input ('Hit enter to run the CAMD demo: ') ;
try
    camd_demo
catch
    disp (lasterr) ;
    fprintf ('CAMD demo failed\n' )
end

input ('Hit enter to run the COLAMD demo: ') ;
try
    colamd_demo
catch
    disp (lasterr) ;
    fprintf ('COLAMD demo failed\n' )
end

input ('Hit enter to run the CCOLAMD demo: ') ;
try
    ccolamd_demo
catch
    disp (lasterr) ;
    fprintf ('CCOLAMD demo failed\n' )
end

fprintf ('\n\n---- CHOLMOD_ACM_TOMS demos complete\n') ;
