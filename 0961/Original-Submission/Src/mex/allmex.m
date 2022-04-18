function allmex( slicot_src )
%
% Function for compiling the SLICOT and SHHEIG sources using the built-in
% FORTRAN compiler in MATLAB and generating and testing the associated
% gateway functions for skew-Hamiltonian/Hamiltonian eigensolvers.
%
% Arguments:
%
% slicot_src - path to the location of the slicot source files, the default
%              path will be used if not provided.
%
% Contributor:
% M. Voigt, Jul. 2014.
%
% Revisions:
% M. Voigt, Jun. 2015.
%
if nargin == 1
    slicotmex( slicot_src );
else
    slicotmex;
end
shheigmex;
makemex;
testmex;