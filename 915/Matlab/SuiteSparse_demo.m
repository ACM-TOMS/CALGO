function SuiteSparse_demo (matrixpath)
%SUITESPARSE_DEMO a demo of all packages in SuiteSparse (SuiteSparseQR subset)
%
% Example:
%   SuiteSparse_demo
%
% See also cholmod, amd, camd, colamd, ccolamd, spqr

% Copyright (c) Timothy A. Davis, Univ. of Florida

if (nargin < 1)
    try
	% older versions of MATLAB do not have an input argument to mfilename
	p = mfilename ('fullpath') ;
	t = strfind (p, filesep) ;
	matrixpath = [ p(1:t(end)) 'CXSparse/Matrix' ] ;
    catch                                                                   %#ok
	% mfilename failed, assume we're in the SuiteSparse directory
	matrixpath = 'CXSparse/Matrix' ;
    end
end

input ('Hit enter to run the CHOLMOD demo: ') ;
try
    cholmod_demo
catch                                                                       %#ok
    disp (lasterr) ;                                                        %#ok
    fprintf ('CHOLMOD demo failed\n' )
end

input ('Hit enter to run the CHOLMOD graph partitioning demo: ') ;
try
    graph_demo
catch                                                                       %#ok
    disp (lasterr) ;                                                        %#ok
    fprintf ('graph_demo failed, probably because METIS not installed\n') ;
end

input ('Hit enter to run the AMD demo: ') ;
try
    amd_demo
catch                                                                       %#ok
    disp (lasterr) ;                                                        %#ok
    fprintf ('AMD demo failed\n' )
end

input ('Hit enter to run the CAMD demo: ') ;
try
    camd_demo
catch                                                                       %#ok
    disp (lasterr) ;                                                        %#ok
    fprintf ('CAMD demo failed\n' )
end

input ('Hit enter to run the COLAMD demo: ') ;
try
    colamd_demo
catch                                                                       %#ok
    disp (lasterr) ;                                                        %#ok
    fprintf ('COLAMD demo failed\n' )
end

input ('Hit enter to run the CCOLAMD demo: ') ;
try
    ccolamd_demo
catch                                                                       %#ok
    disp (lasterr) ;                                                        %#ok
    fprintf ('CCOLAMD demo failed\n' )
end

input ('Hit enter to run the SuiteSparseQR demo: ') ;
try
    spqr_demo
catch                                                                       %#ok
    disp (lasterr) ;                                                        %#ok
    fprintf ('SuiteSparseQR demo failed\n' )
end

fprintf ('\n\n---- SuiteSparse demos complete\n') ;
