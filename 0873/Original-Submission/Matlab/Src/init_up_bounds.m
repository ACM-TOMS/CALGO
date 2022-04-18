%
%   File:              init_up_bounds.m
%
%   Description:       Computes initial alphaU, deltaU
%
%   Authors:           Marielba Rojas
%                        mr@imm.dtu.dk
%                      Sandra Santos
%                      Danny Sorensen
%
%   Version:           1.2
%
%   System:            MATLAB 6.0 or higher 
%
%   Date:              15 March 2007
%
%   Functions called:  diagonal, disp, matvec, min, norm, rand, sprintf
%
%   Call: [alphaU,deltaU] = init_up_bounds(Balpha,norg,Delta,deltaU,message_level);
%

function [alphaU,deltaU] = init_up_bounds(Balpha,norg,Delta,deltaU,message_level)

   n = Balpha.dim - 1;

   if (isa(deltaU,'char'))
      if (strcmp(deltaU,'rayleigh'))
         u           = rand(n,1);
         u           = u/norm(u);
         Balpha.bord = 0;
         deltaU      = u'*matvec(u,Balpha); 
      else
         deltaU = full(min(diag(Balpha.H)));
      end
   end

%
%    Initializes alphaU
%
     alphaU = deltaU + norg*Delta;

     if (message_level == 2)
        disp('--------------')
        disp('init_up_bounds')
        disp('--------------')
        disp(sprintf('\nalphaU = %6.4e',alphaU))
        disp(sprintf('deltaU = %6.4e\n',deltaU))
     end
