function klupackage_install (do_demo)
%klupackage_install: compiles and installs KLU and the packages it requires
% (AMD, BTF, COLAMD).  You must be in the KLUpackage directory to run this file.
%
% Packages in KLUpackage:
%
% AMD            sparse symmetric approximate minimum degree ordering
% COLAMD         sparse column approximate minimum degree ordering
% KLU            sparse LU factorization (left-looking)
% BTF            permutation to block triangular form (like dmperm)
%
% Example:
%    klupackage_install
%    help klupackage_install      % for more details
%
% See also AMD, COLAMD, KLU, BTF, PATHTOOL, PATH.

% KLU and BTF: Copyright 2010, Timothy A. Davis and Eka Palamadai Natarajan.
% COLAMD is co-authored with Stefan Larimore, Esmond Ng, and John Gilbert.
% AMD is co-authored with Patrick Amestoy and Iain Duff.
% http://www.cise.ufl.edu/research/sparse

paths = { } ;
klupackage = pwd ;

% print the introduction
help klupackage_install

% add klupackage to the path
fprintf ('\nPlease wait while KLU, BTF, AMD, and COLAMD are compiled and installed...\n') ;
paths = add_to_path (paths, klupackage) ;

% compile and install AMD
try
    cd ([klupackage '/AMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    amd_make
catch                                                                       %#ok
    disp (lasterr) ;
    fprintf ('AMD not installed\n') ;
end

% compile and install COLAMD
try
    cd ([klupackage '/COLAMD/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    colamd_make
catch                                                                       %#ok
    disp (lasterr) ;
    fprintf ('COLAMD not installed\n') ;
end

% compile and install BTF
try
    cd ([klupackage '/BTF/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    btf_make
catch                                                                       %#ok
    disp (lasterr) ;
    fprintf ('BTF not installed\n') ;
end

% compile and install KLU
have_metis = 0 ;
try
    cd ([klupackage '/KLU/MATLAB']) ;
    paths = add_to_path (paths, pwd) ;
    klu_make (have_metis) ;
catch                                                                       %#ok
    disp (lasterr) ;
    fprintf ('KLU not installed\n') ;
end


% post-install wrapup

cd (klupackage)
fprintf ('KLU, BTF, AMD, COLAMD are now installed.\n') ;

if (nargin < 1)
    % ask if demo should be run
    y = input ('Hit enter to run the KLU demo (or "n" to quit): ', 's') ;
    if (isempty (y))
        y = 'y' ;
    end
    do_demo = (y (1) ~= 'n') ;
end
if (do_demo)
    try
	btf_demo ;
        figure (2) ;
        klu_demo ;
    catch                                                                   %#ok
        disp (lasterr) ;
	fprintf ('KLUpackage demo failed\n') ;
    end
end

fprintf ('\nKLU, BTF, AMD, COLAMD installation is complete.  The following paths\n') ;
fprintf ('have been added for this session.  Use pathtool to add them\n') ;
fprintf ('permanently.  If you cannot save the new path because of file\n');
fprintf ('permissions, then add these commands to your startup.m file.\n') ;
fprintf ('Type "doc startup" and "doc pathtool" for more information.\n\n') ;
for k = 1:length (paths)
    fprintf ('addpath %s\n', paths {k}) ;
end
cd (klupackage)

%-------------------------------------------------------------------------------
function paths = add_to_path (paths, newpath)
% add a path
addpath (newpath) ;
paths = [paths { newpath } ] ;						    %#ok
