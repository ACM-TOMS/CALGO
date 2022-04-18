clc
clear

n=8;

a=0.70005;
b=0.79995;
numPuntos=400;
puntos=linspace(a,b,numPuntos);

x=sym('x','real');
p=sym((x-3/4)^7*(x-1));

disp(' ');
disp('This test file test1.m shows the performance of the Polynomial_1.0 package for the evaluation');
disp(' ');
disp('POLYNOMIAL considered: p(x)=(x-3/4)^7*(x-1)');
disp(' ');
disp('Evaluating p(x)=(x-3/4)^7*(x-1) by Horner algorithm at 400 points');
disp('equally distributed between 0.70005 and 0.79995, both included');
disp('(results are not shown in the screen because of its size)');


coefMon=double(coeffs(p,x));
evalHorner=Horner(coefMon,n,puntos);

disp(' ');
disp('Constructing the polynomial p(x) from its representation in the power'); 
disp('basis with the corresponding constructor of the class (Polynomial(coef,''m''))');
disp('(let us remind that a change of basis is performed in this case)');
disp(' ');

poly = Polynomial(coefMon','m');

disp(' ');
disp('Evaluating p(x) with Eval method of the class Polynomial with the');
disp('default precision 1e-12 (results are not shown in the screen because');
disp('of its size)');
disp(' ');

evalAcc = poly.Eval(puntos);

disp(' ');
disp('Computing the exact values of the polynomial p(x) by using symbolic');
disp('capacities of Matlab (results are not shown in the screen because');
disp('of its size)');
disp(' ');
evalExacta=zeros(1,numPuntos);
for i=1:numPuntos
    evalExacta(i)=double(subs(p,sym(puntos(i))));
end

disp(' ');
disp('Computing the errors for Horner algorithm and Eval method');
ErrorAcc=abs((evalAcc(1,:)-evalExacta)./evalExacta);
ErrorHorner=abs((evalHorner-evalExacta)./evalExacta);
disp(' ');

disp('Mean of the relative errors for the Horner algorithm');
disp(mean(ErrorHorner));
disp('Mean of the relative errors for the Eval method of the class');
disp(mean(ErrorAcc));

disp(' ')
disp('We can observe that Eval method is more accurate than Horner algorithm');
disp('in spite of the change of basis');
disp(' ');
disp(' ');
disp('Press any key to see a graphic with the errors');
pause

figure;
hold off;
semilogy(puntos,ErrorHorner);
hold on;
semilogy(puntos,ErrorAcc,'--');
legend('Horner algorithm','Eval method','Location','SouthEast');
xlabel(texlabel('p(x)=(x-3/4)^7*(x-1)'));
ylabel('Relative errors');
hold off;

