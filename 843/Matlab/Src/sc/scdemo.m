function scdemo
%SCDEMO  Demonstrate the Schwarz-Christoffel Toolbox.

%   Copyright 1997 by Toby Driscoll. 
%   $Id: scdemo.m,v 1.5 1999/01/28 22:17:42 tad Exp $

t = {'Schwarz-Christoffel Toolbox' 'demonstrations'}; % menu title
demo = { 'demtut','deminf','demlong','demmult','demfaber' }; % demos
% Loop until user selects 'quit'
while 1
  k = menu(t,'Tutorial','Infinite vertices','Elongated polygons', ...
      'Faber polynomials','Quit'); 
  if k > 4
    break
  end

  progname = {'scdtutor','scdinf','scdlong','scdfaber'};
  eval(progname{k})
  pause(5)
  h = findall(0,'name','Slideshow Player');
  waitfor(h)
end
