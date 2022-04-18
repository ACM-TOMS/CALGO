function[Psi] =  CalculatePsiConstantsGC(theta,n)
% [Psi] =  CalculatePsiConstantsGC(theta,n)
% CalculatePsiConstantsGC:  Calculates the total Psi normalization constant from the Theta Parameter and psi-s 
% vector for the Cayley distance.
% 
% INPUTS
% theta: The theta parameter
% n:  Population size 
% OUTPUTS
% Psi: Calculated psi values.
%
% Created version 04/15/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/15/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

j = [1:n-1]; % Index from 1..n-1 according to the papper
Psi = (n-j) .* exp(-theta) + 1;


