function[invperm] = Invert(perm)
% [invperm] = Invert(perm)
% Invert:  Inverts a permutation.
%             
% perm:  The permutation to invert
% OUTPUTS
% invperm: The inverted permutation
%
% Created version 12/22/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

[auxval,invperm] =sort(perm);