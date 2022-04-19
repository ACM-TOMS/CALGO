function [ VV ] = tjsr_getpolytope( type );   
% [ VV ] = tjsr_getpolytope( type );   
% This function belongs to tjsr!
% select the vertices which make up the current polytope
% This is not a polytope which is suitable to compute a norm with (in general)
%
% Input:
%   type            the type struct from tjsr
%
% Output:
%   VV              selected vertices
%
% Written by tommsch, 2018

    VV = [type.cyclictree.V{:}]; %get all vertices
    idx = true(1,size(VV,2)); %indices of all vertices
    
    idx_norm = [type.cyclictree.norm{:}]>1-type.opt.epspolytope; 
    idx = idx & idx_norm; %only select vertices which are outside or at the border of the polytope.-
    
    val = min(type.opt.simplepolytope, (type.JSR(2)/type.JSR(1)-1)/1000)-type.opt.epslinprog;
    idx_dist = [type.cyclictree.norm{:}] > 1+val;
    idx = idx & idx_dist; %only select vertices which are far away from other vertices 
    
    VV = VV(:,idx); %remove all unselected vertices
    VV = removezero(VV,2); %remove empty columns
    
end