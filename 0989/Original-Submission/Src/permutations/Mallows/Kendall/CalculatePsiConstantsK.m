function[Psi] =  CalculatePsiConstantsK(theta,n)
% [Psi] =  CalculatePsiConstantsK(theta,n)
% CalculatePsiConstantsK:  Calculates the total Psi normalization constant from the Theta Parameter
% and psi-s vector for the Kendall distance.
%         
% INPUTS
% theta: The theta parameter
% n:  Population size 
% OUTPUTS
% Psi: Calculated psi values.
%
% Created version 12/11/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

j = [1:n-1]; % Index from 1..n-1 according to the papper
Psi = (1 - exp((-1)*(n-j+1)*(theta))) ./ (1 - exp((-1)*theta));


