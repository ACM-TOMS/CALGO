function makemex
%
% Function for generating gateway functions from MEX files for 
% skew-Hamiltonian/Hamiltonian eigensolvers.
%
% Contributor:
% M. Voigt, Jul. 2013.
%
% Revisions:
% M. Voigt, Jul. 2014, Jun. 2015.
%
%% Set flags depending on machine architecture.
%
flags = '';
is64 = ( ~isempty( strfind( computer, '64' ) ) );
if ( is64 )
    %
    % 64-bit MATLAB
    %
    flags = 'FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer -fdefault-integer-8" -largeArrayDims';
end
%
%% Set location of the SLICOT and SHHEIG libraries. You might want to modify these. 
%
libslicot = '../mexslicot.a';
libshheig = '../mexshheig.a';
%
%% Mexing.
%
shheig_mex = {
    'skewHamil2eig', ...
    'skewHamil2eig_nb', ... 
    'skewHamil2feig', ...
    'skewHamildefl', ...
    'skewHamildeflf', ...
    'skewHamildeflfZ', ... 
    'skewHamildefl_nb', ... 
    'skewHamildeflZ', ...
    'skewHamildeflZ_nb', ... 
    'skewHamileig', ...
    'skewHamileig_nb', ... 
    'symplURV', ...
    };
%
for k = 1:length(shheig_mex)
    file = shheig_mex{k};
    fprintf( 'mex %s %s.F %s %s -lmwlapack -lmwblas\n', flags, file, libshheig, libslicot );
    eval( sprintf( 'mex %s %s.F %s %s -lmwlapack -lmwblas\n', flags,  file, libshheig, libslicot ) );
end
