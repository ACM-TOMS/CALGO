function[Psi] =  CalculatePsiConstantsGK(theta,n)
% [Psi] =  CalculatePsiConstantsGK(theta,n)
% CalculatePsiConstantsGK:  Calculates the total Psi normalization constant from the Theta vector 
% parameter and psi-s vector for the Kendall distance.
%         
% INPUTS
% theta: The theta vector parameter
% n:  Population size 
% OUTPUTS
% Psi: Calculated psi values.
%
% Created version 04/11/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/15/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

j = [1:n-1]; % Index from 1..n-1 according to the papper
Psi = (1 - exp((-1)*(n-j+1).*theta)) ./ (1 - exp((-1).*theta));


