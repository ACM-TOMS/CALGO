% Example 1: Computing the nodes and weight in an n-point rational Fejer 
%            quadrature formula based on a sequence of (n-1) poles

clear; 
close; 
clc; 

disp('The sequence of poles is [2:1:10].');
sgl = [2:1:10];
[x,err] = rfejer(sgl);


disp(sprintf('\n The quadrature formula should be exact for integrals of the form:'));
disp('int_{-1}^{1} ( (x-a).^(-1) ) dx, for a = 2,...,10.');
ErrExact = max(err(2,:));
disp(sprintf('Exact maximal relative error on the approximations: %1.16e',ErrExact));


disp(sprintf('\n We now reverse the order of the poles:'));
disp('The new sequence of poles is [10:-1:2].');
sgl2 = sgl(end:-1:1);
[y,err2] = rfejer(sgl2);


disp(sprintf('\n The quadrature formula should be exact for the same integrals as before.'));
ErrExact2 = max(err2(2,:));
disp(sprintf('Exact maximal relative error on the approximations: %1.16e',ErrExact2));


disp(sprintf('\nTheoretically the nodes and weights do not depend on the')); 
disp(sprintf('order of the poles. Hence, the weights should be identical'));
disp(sprintf('for both sequences of poles.'));
err3 = abs(x(2,:)-y(2,:))./abs(x(2,:));
ErrExact3 = min(err3);
disp(sprintf('\n In practice, however, the weights can contain large errors.'));
disp(sprintf('When comparing the computed weights for both sequences of poles,'));
disp(sprintf('we obtain a minimal relative distance: %1.16e',ErrExact3));


disp(sprintf('\n Despite the large errors, both quadrature formulae perform'));
disp(sprintf('equally well for the approximation of integrals of the form:'));
disp('int_{-1}^{1} f(x) dx, where the function f is arbitrary.');
fun = @(x) 1./(x-1.8);
Q1 = fun(x(1,:))*x(2,:)';
Q2 = fun(y(1,:))*y(2,:)';
ExactValue = log(1-1.8)-log(-1-1.8);
E1 = abs(Q1-ExactValue)/abs(ExactValue);
E2 = abs(Q2-ExactValue)/abs(ExactValue);
E3 = abs(Q1-Q2)/abs(Q1);
disp(sprintf('\n Consider the case in which f(x) = (x-1.8).^(-1).'));
disp(sprintf('Note that none of poles coincides with the pole of f at x = 1.8.'))
disp(sprintf('The exact value for the integral: %1.16e',ExactValue));
disp(sprintf('The approximation obtained from the first sequence of poles: %1.16e',Q1));
disp(sprintf('with exact relative error: %1.16e',E1));
disp(sprintf('The approximation obtained from the second sequence of poles: %1.16e',Q2));
disp(sprintf('with exact relative error: %1.16e',E2));
disp(sprintf('The relative distance between the two approximations: %1.16e',E3));
