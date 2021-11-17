%   Pc_sorted = cmpoles(P, Pc);
%
%   Sorts the vector Pc to be in the order of P.  
%
%   This function is useful to compare, e.g., 
%   the desired and assigned poles of a system.  

%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.0, Copyright (c) 2002-2009 NICONET e.V.
%
%   Based on function localcmpeig, included in place.m 
%   31 July, 2002

function Pc_sorted = cmpoles(P, Pc)

Pc_sorted = zeros( size( P ) );

for i = 1 : length( P ),
    [dummy, j] = min( abs( P(i) - Pc(:) ) );
    Pc_sorted(i,1) = Pc(j);  Pc(j) = [];
end
