function[fthetadev] =  UlamThetaDevFunction(theta,newton_params)
% Ulam Theta Dev parameter estimation function.
%
% Created version 02/12/2015. Roberto Santana (roberto.santana@ehu.es) 
%
% Last version 04/20/2015. Roberto Santana (roberto.santana@ehu.es)  
%

n = newton_params{1};              %number of variables
ncount = newton_params{2};

vectornjs = [1:n-1]; % 1..n-1
aux  = ncount(1:n-1) ./ exp(theta*vectornjs); 
aux_numer1 = (aux .* vectornjs)  .* vectornjs;
aux_numer3 = (aux .* vectornjs);

numer1 = sum(aux_numer1);
numer2 = sum(aux);
numer3 = sum(aux_numer3);

fthetadev = (-numer1*numer2 - numer3*numer3)/(numer2*numer2);
