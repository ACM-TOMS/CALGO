%
%   File:              upd_deltaU.m
%
%   Description:       updates approximation to smallest eigenpair of H
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
%   Functions called:  disp, min, sprintf
%
%   Call: deltaU = upd_deltaU(epair,g,deltaU,message_level);
%

function [deltaU] = upd_deltaU(epair,g,deltaU,message_level)

     lambda = epair.lambda; 
     nu     = epair.nu; 
     u      = epair.u; 

     deltaN     = lambda - nu * (g'*u)/(1 - nu^2);
     old_deltaU = deltaU;
     deltaU     = min([deltaU,deltaN]);
  
     if (message_level == 2)
        disp('----------')
        disp('upd_deltaU')
        disp('----------')
        if (old_deltaU ~= deltaU)
           disp(sprintf('\nupdating deltaU'));
        else
           disp(sprintf('\nnot updating deltaU'));
        end
        disp(sprintf('\n   deltaN        deltaU'))
        disp(sprintf('%e %e\n',deltaN,old_deltaU))
     end
