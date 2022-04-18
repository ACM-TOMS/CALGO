%
%   File:              safe_alpha1.m
%
%   Description:       Safeguards alpha1
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
%   Call: alpha = safe_alpha1(alpha,alphaL,alphaU,deltaU,ipoint,message_level);
%

function [alpha] = safe_alpha1(alpha,alphaL,alphaU,deltaU,ipoint,message_level)

    alpha1 = alpha;

    lambda = ipoint.lambda; 
    fi     = ipoint.fi;
    norx   = ipoint.norx;

    if ((alpha < alphaL) | (alpha > alphaU))
       alpha = deltaU + fi + (norx^2)*(deltaU - lambda);
       if (alpha < alphaL | alpha > alphaU)
          alpha2 = alpha;
          alpha  = (alphaL + alphaU)/2;
       end
    end

    if (message_level == 2)
       disp('-----------')
       disp('safe_alpha1')
       disp('-----------')
       if (alpha ~= alpha1)
          disp(sprintf('\nsafeguarding alpha1'))
	  if (exist('alpha2'))
             if (alpha ~= alpha2)
                disp(sprintf('using bisection!!'))
             end 
          end 
       else
          disp(sprintf('\nnot safeguarding alpha1'))
       end
       disp(sprintf('alpha = %e\n',alpha))
    end
