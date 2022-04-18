% Example script to test AMGKQ in two dimensions.
% Copyright (C) 2014 Robert W. Johnson
% Based on TEST_INT_2D (C) 2009 & 2011 John Burkardt
%
% This file is part of AMGKQ.  See AMGKQ.M for details.
% Copyright (C) 2014 Robert W. Johnson, Alphawave Research.
% This is free software; see GPLv3 or greater for copying conditions.
% There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  For details, see GPLv3 or greater.

% preparations
close all
clear all
tabcell = cell(1,9);
ntab = 0;
disp(' ')

fx = @(x) 1 ./ ( 1 - x(1,:) .* x(2,:) )
a(1:2) = 0
b(1:2) = 1
exact = pi^2 / 6
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) 1 ./ sqrt ( 1 - x(1,:).^2 .* x(2,:).^2 )
a(1:2) = -1
b(1:2) =  1
exact = 2 * pi * log ( 2 )
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) 1 ./ sqrt ( 2 - x(1,:) - x(2,:) )
a(1:2) = -1
b(1:2) =  1
exact = ( 16 / 3 ) * ( 2 - sqrt ( 2 ) )
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) 1 ./ sqrt ( 3 - x(1,:) - 2 * x(2,:) )
a(1:2) = -1
b(1:2) =  1
exact = ( sqrt ( 32 ) / 3 ) * ( sqrt ( 27 ) - sqrt ( 8 ) - 1 )
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) sqrt ( x(1,:) .* x(2,:) )
a(1:2) = 0
b(1:2) = 1
exact = 4 / 9
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) exp ( - ( ( x(1,:) - 4 ).^2 + ( x(2,:) - 1 ).^2 ) )
a(1:2) = 0
b(1:2) = 5
exact = 0.25 * pi * ( erf ( 1 ) + erf ( 4 ) ).^2
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) abs ( x(1,:).^2 + x(2,:).^2 - 0.25 )
a(1:2) = -1
b(1:2) =  1
exact = ( 5 / 3 ) + ( pi / 16 )
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) sqrt ( abs ( x(1,:) - x(2,:) ) )
a(1:2) = 0
b(1:2) = 1
exact = 8 / 15
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

% redo difficult integrands
numfx = ntab
redofx = find(abs(cell2mat(tabcell(:,9))) > 1.5e-8)
numredo = length(redofx)
for kfx = redofx'
[fx, a, b, exact] = deal(tabcell{kfx,1:4})
[RES, ERR, NSUB, FL] = amgkq(fx,a',b',[],0,1e3), ACC = RES-exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};
end

% format table for display
fnos = [1:ntab]';
fnos = fnos - floor((fnos - 1)/numfx)*numredo;
headstr = 'No.               A                 B    EXACT      RES      ERR     NSUB       FL      ACC';
tabshow = strvcat(headstr,num2str([fnos,cell2mat(tabcell(:,2:end))],'%8.2g '))

