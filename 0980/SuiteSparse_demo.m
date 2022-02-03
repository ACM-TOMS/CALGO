function SuiteSparse_demo (matrixpath, dopause)
%SUITESPARSE_DEMO a demo of all packages in SuiteSparse (SPQR subset)
%
% Example:
%   SuiteSparse_demo
%
% See also umfpack, cholmod, amd, camd, colamd, ccolamd, btf, klu, spqr,
%   CSparse, CXSparse, ldlsparse

% Copyright 2016, Timothy A. Davis, http://www.suitesparse.com.

if (nargin < 1 || isempty (matrixpath) || ~ischar (matrixpath))
    try
	% older versions of MATLAB do not have an input argument to mfilename
	p = mfilename ('fullpath') ;
	t = strfind (p, '/') ;
	matrixpath = [ p(1:t(end)) 'CXSparse/Matrix' ] ;
    catch me    %#ok
	% mfilename failed, assume we're in the SuiteSparse directory
	matrixpath = 'CXSparse/Matrix' ;
    end
end

if (nargin < 2)
    dopause = false ;
end

if (dopause)
    input ('Hit enter to run the CHOLMOD demo: ', 's') ;
end
try
    cholmod_demo
catch me
    disp (me.message) ;
    fprintf ('CHOLMOD demo failed\n' )
end

if (dopause)
    input ('Hit enter to run the CHOLMOD graph partitioning demo: ', 's') ;
end
try
    graph_demo
catch me
    disp (me.message) ;
    fprintf ('graph_demo failed, probably because METIS not installed\n') ;
end

if (dopause)
    input ('Hit enter to run the AMD demo: ', 's') ;
end
try
    amd_demo
catch me
    disp (me.message) ;
    fprintf ('AMD demo failed\n' )
end

if (dopause)
    input ('Hit enter to run the CAMD demo: ', 's') ;
end
try
    camd_demo
catch me
    disp (me.message) ;
    fprintf ('CAMD demo failed\n' )
end

if (dopause)
    input ('Hit enter to run the COLAMD demo: ', 's') ;
end
try
    colamd_demo
catch me
    disp (me.message) ;
    fprintf ('COLAMD demo failed\n' )
end

if (dopause)
    input ('Hit enter to run the CCOLAMD demo: ', 's') ;
end
try
    ccolamd_demo
catch me
    disp (me.message) ;
    fprintf ('CCOLAMD demo failed\n' )
end

if (dopause)
    input ('Hit enter to run the SuiteSparseQR demo: ', 's') ;
end
try
    spqr_demo
catch me
    disp (me.message) ;
    fprintf ('SuiteSparseQR demo failed\n' )
end

fprintf ('\n\n---- SuiteSparse demos complete\n') ;
