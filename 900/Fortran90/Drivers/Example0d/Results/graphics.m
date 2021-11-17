% clear and close
aux = sprintf('cleaning and closing ...\n') ; disp( aux );
clear all;
close all;

% load info
aux = sprintf('loading files ...\n') ; disp( aux );
load truth.txt truth
load model.txt model
load ekf.txt ekf
load rrsqrtkf.txt rrsqrtkf
load enkf.txt enkf
load rrsqrtenkf.txt rrsqrtenkf
load obser.txt obser

% truth
aux = sprintf('building true vector ...\n') ; disp( aux );
xt = truth( : , 1 );
yt = truth( : , 2 );

% model
aux = sprintf('building model vector ...\n') ; disp( aux );
xm = model( : , 1 );
ym = model( : , 2 );

% observations
aux = sprintf('building observation vector ...\n') ; disp( aux );
xo = obser( : , 1 );
yo = obser( : , 2 );

% ekf
aux = sprintf('building ekf vector ...\n') ; disp( aux );
xekf = ekf( : , 1 );
yekf = ekf( : , 2 );

% rrsqrtkf
aux = sprintf('building rrsqrtkf vector ...\n') ; disp( aux );
xrrsqrtkf = rrsqrtkf( : , 1 );
yrrsqrtkf = rrsqrtkf( : , 2 );

% enkf
aux = sprintf('building enkf vector ...\n') ; disp( aux );
xenkf = enkf( : , 1 );
yenkf = enkf( : , 2 );

% rrsqrtenkf
aux = sprintf('building rrsqrtenkf vector ...\n') ; disp( aux );
xrrsqrtenkf = rrsqrtenkf( : , 1 );
yrrsqrtenkf = rrsqrtenkf( : , 2 );

% plotting
aux = sprintf('plotting ...\n') ; disp( aux );
figure( 1 );
hold on
plot( xt , yt , 'g-' , 'LineWidth' , 3 );
plot( xm , ym , 'b-' );
plot( xo , yo , 'k*' );
plot( xekf , yekf , 'r-' );
plot( xrrsqrtkf , yrrsqrtkf , 'b-.' )
plot( xenkf , yenkf , 'c-' , 'color' , [0.645 0.164 0.164] );
plot( xrrsqrtenkf , yrrsqrtenkf , 'm-' );
hold off
h = get( gcf , 'children' );
set( h , 'fontsize' , 10 );
xlabel( 'TIME' );
ylabel( ' VALUES OF' );
legend( 'TRUTH' , 'MODEL' , 'OBSERVATIONS' , 'EKF' , 'RRSQRTKF' , 'ENKF' , 'RRSQRTENKF' );
grid on

% plotting
aux = sprintf('saving figure in example0d_values.eps ...\n') ; disp( aux );
print -depsc example0d_values.eps


