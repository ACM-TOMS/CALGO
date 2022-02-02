function spdemo
% SPDEMO   A 2D-example for multi-linear sparse grid
% interpolation using the Clenshaw-Curtis grid and vectorized
% processing of the model function.
%
%    See also SPINTERP, SPVALS.

% Author : Andreas Klimke, Universit�t Stuttgart
% Version: 1.1
% Date   : September 29, 2003

% Some function f
f = inline('1./((x*2-0.3).^4 +(y*3-0.7).^2+1)');

% Define problem dimension
d = 2;

% Create full grid for plotting
gs = 33;
[X,Y] = meshgrid(linspace(0,2,gs),linspace(-1,1,gs));

% Set options: Switch vectorized processing on.
options = spset('Vectorized', 'on');

% Compute sparse grid weights over domain [0,2]x[-1,1]
z = spvals(f, d, [0 2; -1 1], options);

% Compute inpterpolated values at full grid
ip = spinterp(z, X, Y);

% Plot original function, interpolation, and error
subplot(1,3,1);
mesh(X,Y,f(X,Y));
title('original');

subplot(1,3,2);
mesh(X,Y,ip);
title('interpolated');

subplot(1,3,3);
mesh(X,Y,abs(f(X,Y)-ip));
title('absolute error');

disp(' ');
disp('Sparse grid representation of the function:');
z