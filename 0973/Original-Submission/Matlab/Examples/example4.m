% Example 4: Numerical approximation of an integral on the interval [0,1]  
%            where the integrand is a vector-valued function with different
%            essential singularities

clear; 
close; 
clc; 

disp(sprintf('Vector-valued function with distinct singularities.'));
disp('int_{0}^{1} sin(1/(x-k*omeg )) dx for k=1:.1:2 and omeg=1.1');

% The integrand does not accept a vector argument when k is not a scalar 
fun = @(x,w,k) sin(1./(x-k*w));
omeg = 1.1;

disp(sprintf('\nThe integrand has essential singularities'));
disp(sprintf('at k* %1.3f for k=1, 1.1, 1.2, ... ,2.',omeg));
S = [1:.1:2];
fx = @(x) fun(x,omeg,S);

sizeFx = size(fx(0));
Exact = zeros(sizeFx);
Exact(:) = [-.594506278435972795476222305809, -.680954390598676563313531763785, ...
    -.821216684622060163075117363820, -.847197587427479419618061638019, ...
    -.820172814035944885000378601742, -.775892308741661865885845669493, ...
    -.728238071367740968289963827901, -.682324585956527728754834805052, ...
    -.639873232845142628666014318056, -.601271625610329305498003832839, ...
    -.566388120017313046000120829223];


disp(sprintf('\nFirst we consider poles that coincide with the singularities.'));
sgl1 = omeg*[S,S,S];
n = length(sgl1)+1;
disp(sprintf('The maximal number of iterations is %1.3f',n));

% First we need to map the interval [0,1] onto the interval [-1,1]
[fout,sglout ] = transf( fx , sgl1 , [0,1] );

% Since the integrand does not accept a vector argument, we need to use the 
% third input argument of rfejer and set it to TRUE  
[NumInt1,Err1] = rfejer(sgl1,fout,'Array',true);
ErrMax = max(Err1(:));
disp(sprintf('Estimated maximal relative error on the computed values: %1.16e',ErrMax));
ErrExact = abs(Exact-NumInt1)./abs(Exact);
ErrExact = max(ErrExact(:));
disp(sprintf('Exact maximal relative error : %1.16e',ErrExact));

omeg2 = omeg*mean(S);
disp(sprintf('\nNext we consider a multiple pole at %1.3f.',omeg2));
sgl2 = omeg2*ones(1,n-1);
disp(sprintf('The maximal number of iterations is again %1.3f',n));

% First we need to map the interval [0,1] onto the interval [-1,1]
[fout,sglout ] = transf( fx , sgl2 , [0,1] );

% Since the integrand does not accept a vector argument, we need to use the 
% third input argument of rfejer and set it to TRUE  
[NumInt2,Err2] = rfejer(sgl2,fout,'Array',true);
ErrMax = max(Err2(:));
disp(sprintf('Estimated maximal relative error on the computed values: %1.16e',ErrMax));
ErrExact = abs(Exact-NumInt2)./abs(Exact);
ErrExact = max(ErrExact(:));
disp(sprintf('Exact maximal relative error : %1.16e',ErrExact));
