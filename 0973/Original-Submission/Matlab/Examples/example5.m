% Example 5: Numerical approximation of an integral on the interval [-1,1]  
%            for which the singularities of the integrand are not known
%            explicitly

clear; 
close; 
clc; 

disp(sprintf('Integrand with unknown singularities.'));
disp('int_{-1}^{1} exp(-x^2) dx');
disp('The integrand tends to infinity when the imaginary part of x tends');
disp('to infinity. Hence, a rational quadrature rule with all poles at');
disp('infinity (i.e., a classical polynomial quadrature rule) will work');
disp('perfectly well for this example. But let us assume for a moment that');
disp('the exact location of the singularities are unknown. One usually');
disp('relies then on the singularities of a [n/m] Pade approximant, where');
disp('n and m denote the degree of the numerator and denominator polynomial');
disp('respectively.')

fun = @(x) exp(-x.^2);

Exact = 1.49364826562485405079893487226;

disp('For the poles we use the zeros of the Maclaurin polynomial of');
disp('degree 20 of exp(x^2):');
disp('1+x^2+(1/2)*x^4+(1/6)*x^6+(1/24)*x^8+(1/120)*x^10+(1/720)*x^12+');
disp('    +(1/5040)*x^14+(1/40320)*x^16+(1/362880)*x^18+(1/3628800)*x^20');
disp('(which corresponds to using the singularities of the [1/20] Pade');
disp('approximant of exp(-x^2)).')

sgl = [2.2288448993805862124+1.2620932441870966540*i, ...
-2.2288448993805862124-1.2620932441870966540*i, ... 
1.5865563357512161835+1.5655540425711388648*i, ...
-1.5865563357512161835-1.5655540425711388648*i, ...
1.0810988862502763130+1.7436842639282667696*i, ... 
-1.0810988862502763130-1.7436842639282667696*i, ...
.63185494016061721900+1.8479113712057029298*i, ...
-.63185494016061721900-1.8479113712057029298*i, ...
.20811231677867307823+1.8966250895533320338*i, ...
-.20811231677867307823-1.8966250895533320338*i];

n = length(sgl)+1;
disp(sprintf('The maximal number of iterations is %1.3f',n));
[NumInt,Err] = rfejer(sgl,fun);
disp(sprintf('Computed value: %1.16e',NumInt));
disp(sprintf('Estimated relative error : %1.16e',Err));
ErrExact = abs(Exact-NumInt)/abs(Exact);
disp(sprintf('Exact relative error : %1.16e',ErrExact));
