%
%   File:              inter_point.m
%
%   Description:       Selects interpolation point
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
%   Call: ipoint = inter_point(epair1,epair2,iterate,g,message_level);
%

function ipoint = inter_point(epair1,epair2,iterate,g,message_level)

         if (iterate == 0)
            ipoint = [];
         elseif (iterate == 1)
            x    = epair1.u/epair1.nu;
            norx = epair1.noru/epair1.anu;

            ipoint.lambda = epair1.lambda;
            ipoint.fi     = - g'*x;
            ipoint.norx   = norx;
         elseif (iterate == 2)
            x    = epair2.u/epair2.nu;
            norx = epair2.noru/epair2.anu;

            ipoint.lambda = epair2.lambda;
            ipoint.fi     = - g'*x;
            ipoint.norx   = norx;
         end

         if (message_level == 2)
            disp('-----------');
            disp('inter_point');
            disp('-----------');
            if (iterate == 0)
               disp(sprintf('\nNo iterate!'));
            elseif (iterate == 1)
               disp(sprintf('\nusing 1st'));
            else
               disp(sprintf('\nusing 2nd'));
            end
            disp(' ');
         end
