%
%   File:              upd_alpha_safe.m
%
%   Description:       Updates safeguards for alpha
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
%   Functions called:  disp, sprintf
%
%   Call: [alphaL,alphaU] = ...
%          upd_alpha_safe(alpha,alphaL,alphaU,iterate,ipoint,Delta,message_level);
%

function [alphaL,alphaU] = ...
          upd_alpha_safe(alpha,alphaL,alphaU,iterate,ipoint,Delta,message_level)

      norx = ipoint.norx;

      if (iterate == 1)
         if (norx < Delta)
            alphaL = alpha;
         elseif (norx > Delta)
            alphaU = alpha;
         end
      else
         alphaU = alpha;
      end

      if (message_level == 2)
         disp('--------------')
         disp('upd_alpha_safe')
         disp('--------------')
         disp(sprintf('\nalphaL = %6.4e',alphaL));
         disp(sprintf('alphaU = %6.4e\n',alphaU));
      end
