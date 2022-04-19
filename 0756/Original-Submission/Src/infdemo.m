more off
echo on
clc
% This script demonstrates the interpretation of infinite vertices.

pause     % Strike any key to begin (Ctrl-C to abort)

% For the purposes of the Schwarz-Christoffel Toolbox, a polygon
% is represented by two vectors: w, the vertices in positively 
% oriented order, and beta, the corresponding "turning angles."
% At a finite vertex w(j), the meaning of beta(j) is simple:
% pi*beta(j) is the exterior turning angle of the polygon
% at w(j), with a minus sign for left (counterclockwise) turns.  
% By convention, beta(j) is +1 at a slit, so at a finite vertex,
% -1 < beta(j) <= 1.

% If w(j) is infinite, the formal definition of beta(j) is more 
% cumbersome: pi*beta(j) is -2*pi plus the exterior angle formed by 
% the two sides incident on w(j) as they are extended *away* from 
% infinity.  But the interpretation of pi*beta(j) as "turn" is still
% valid.  It turns out that -3 <= beta(j) <= -1 at an infinite vertex.

% Following are some examples that should help clarify matters.

pause     % Strike any key to continue
echo off

figure(gcf)
cla
axis square
axis([-2 2 -2 2])
hold on

int = fill([2,0,0,2],[1,1,-1,-1],[.6,.6,.6]);
edges = plot([2+i,i,-i,2-i],'y-');
plot([-i,i],'.','marker',12)
t1 = text(-.1,1,'-1/2','hor','right','ver','mid');
t2 = text(-.1,-1,'-1/2','hor','right','ver','mid');
t3 = text(2.1,0,'-1','hor','left','ver','mid');

clc
disp(' ')
disp('Here is a three-vertex polygon.  The values of beta are shown')
disp('next to their associated vertices.  The infinite vertex has a')
disp('turn of -1, which is the least possible.')
disp(' ')
disp('  Strike any key to continue')
pause

set(t1,'string','0');
set(t3,'string','-3/2');
set(edges,'xdata',[0,0,2],'ydata',[2,-1,-1])
set(int,'xdata',[0,0,2,2],'ydata',[2,-1,-1,2])

disp(blanks(2)')
disp('The turn at infinity is now -3/2.  Note that the sum of the')
disp('turns is still -2.')
disp(' ')
disp('  Strike any key to continue')
pause


set(t1,'string','1/2','pos',[-.1,.9],'ver','top')
set(t3,'string','-2')
set(edges,'xdata',[-2,0,0,2],'ydata',[1,1,-1,-1])
set(int,'xdata',[-2,0,0,2,2,-2],'ydata',[1,1,-1,-1,2,2])

disp(blanks(2)')
disp('A turn of -2 isn''t really a turn at all.  But the edge returning')
disp('from infinity doesn''t have to be colinear with the outgoing edge.')
disp(' ')
disp('  Strike any key to continue')
pause


set(t1,'pos',[-.1,.4],'string','3/4')
set(t3,'string','-9/4')
set(edges,'xdata',[-2,0,0,2],'ydata',[-2,1,-1,-1])
set(int,'xdata',[-2,0,0,2,2,-2],'ydata',[-2,1,-1,-1,2,2])

disp(blanks(2)')
disp('Turning past -2 produces a more "open" region.')
disp(' ')
disp('  Strike any key to continue')
pause


set(t1,'string','1/2','pos',[.1,.9],'hor','left')
set(t2,'string','1/2','pos',[.1,-.9],'hor','left','ver','bot')
set(t3,'string','-3')
set(edges,'xdata',[2,0,0,2],'ydata',[1,1,-1,-1])
set(int,'xdata',[2,0,0,2,2,-2,-2,2],'ydata',[-1,-1,1,1,2,2,-2,-2])

disp(blanks(2)')
disp('Finally, a turn of -3 is the most allowed.  The interior of this')
disp('polygon is the complement of the first example, with a turn of -1.')
disp(' ')
disp('  Strike any key to continue')
pause

disp(' ')
disp('Here is what a disk map of this region looks like....')
disp(' ')
echo on
w = [-i; i; Inf];
beta = [.5; .5; -3];
[z,c] = dparam(w,beta);
[z,c] = dfixwc(w,beta,z,c,-1+i);
cla
dplot(w,beta,z,c)

echo off     % End of demo

