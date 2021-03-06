function shheigmex
%
% Function for compiling SHHEIG subroutines using the built-in FORTRAN
% compiler of MATLAB. The archive file will be created inte parent
% directory.
%
% Contributor:
% M. Voigt, Jul. 2014.
%
% Revisions:
% -
%
%% Set flags depending on machine architecture.
%
flags = '';
is64 = ( ~isempty( strfind( computer, '64' ) ) );
if ( is64 )
    % 64-bit MATLAB
    flags = 'FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer -fdefault-integer-8" -largeArrayDims';
end
%
%% Set location of the SHHEIG source files.
%
shheig_src = '../src/';
%
%% Compiling.
%
shheig = {
    'DGHFDF', ...
    'DGHFET', ...
    'DGHFEX', ...
    'DGHFEY', ...
    'DGHFST', ...
    'DGHFXC', ...
    'DGHFYR', ...
    'DGHUDF', ...
    'DGHUDP', ...
    'DGHUET', ...
    'DGHUEX', ...
    'DGHUEY', ...
    'DGHURV', ...
    'DGHUSP', ...
    'DGHUST', ...
    'DGHUTP', ...
    'DGHUTR', ...
    'DGHUXC', ...
    'DGHUXP', ...
    'DGHUYR', ...
    'DLACPV', ...
    'ZGHFDF', ...
    'ZGHFEX', ...
    'ZGHFEY', ...
    'ZGHFXC', ...
    'ZGHUDF', ...
    'ZGHUDP', ...
    'ZGHUEX', ...
    'ZGHUEY', ...
    'ZGHUXC', ...
    'ZGHUXP', ...
    'ZLACPV', ...
    };
obj = [];
for k = 1:length(shheig)
    file = shheig{k};
    fprintf( 'mex -c %s %s%s.f\n', flags, shheig_src, file );
    eval( sprintf( 'mex -c %s %s%s.f\n', flags, shheig_src, file ) );
    obj = [ obj ' ' file '.o' ];
end
%
% Archive the object files.
%
unix( 'ar r ../mexshheig.a *.o' );
unix( 'rm *.o' );
