clc
clear

matlab = ver('MATLAB');
matlabVersion = regexp(matlab.Version,'\.','split');
majorVersion = str2num(matlabVersion{1});
minorVersion = str2num(matlabVersion{2});
if majorVersion<7 || (majorVersion==7 && minorVersion<=1)
    error('Version of Matlab must be 7.2 at least');
end

n=21;

a=0.40005;
b=0.59995;
numPuntos=400;
puntos=linspace(a,b,numPuntos);

x=sym('x','real');
r=sym(x*(x-1/2)^20);

disp(' ');
disp('This test file shows the performance of the Polynomial_1.0 package for the evaluation');
disp(' ');
disp('POLYNOMIAL considered: r(x)=x*(x-1/2)^20');
disp(' ')
disp(' ')
disp(' ');
disp('Computing the exact coefficients of r(x)=x*(x-1/2)^20 with respect');
disp('to the Bernstein basis');

aux=sym(zeros(1,n+1));
for i=1:n+1
    aux(i)=sym(nchoosek(n,i-1))*x^(i-1)*(1-x)^(n+1-i);
end

c = sym('c', [1 n+1]);
if (majorVersion==7 && minorVersion >= 14) || majorVersion>=8
    exp=taylor(expand(sum(c.*aux)),'Order',n+1);
    rhs=coeffs(taylor(expand(r),'Order',n+1),x);
else
    exp=taylor(expand(sum(c.*aux)),n+1);
    rhs=coeffs(taylor(expand(r),n+1),x);
end

rhs=[0,rhs];
lhs=coeffs(exp,x);
ecuacion=rhs-lhs;
ExactCoeffBer=solve(ecuacion);

for i=1:n+1
    coeffBer(i)=double(eval(sprintf('ExactCoeffBer.c%d',i)));
    coeffVS(i)=nchoosek(n,i-1)*coeffBer(i);
end



disp(' ');
disp('Evaluating r(x)=x*(x-1/2)^20 by VS algorithm at 400 points equally');
disp('distributed between 0.40005 and 0.59995, both included');
disp('(results are not shown in the screen because of its size)');

[y_VS,errBound_VS] = Vs(coeffVS,puntos);

disp(' ');
disp('Evaluating r(x)=x*(x-1/2)^20 by Casteljau algorithm at 400 points');
disp('equally distributed between 0.40005 and 0.59995, both included');
disp('(results are not shown in the screen because of its size)');

[y_Cast,errBound_Cast] = Casteljau(coeffBer,puntos);

disp(' ');
disp('Constructing the polynomial r(x)=x*(x-1/2)^20 from its exact coefficients'); 
disp('respect to the Bernstein basis with the corresponding constructor of the');
disp('class (Polynomial(coef,''c''))');

poly = Polynomial(coeffBer,'c');

disp(' ');
disp('Evaluating r(x)=x*(x-1/2)^20 by Eval method at 400 points equally');
disp('distributed between 0.40005 and 0.59995 with e pretended precision');
disp('of the epsilon of the machine');
disp('(results are not shown in the screen because of its size)');

y_Eval = poly.Eval(puntos,eps);

disp(' ');
disp('Computing the exact values of the polynomial r(x) by using symbolic');
disp('capacities of Matlab (results are not shown in the screen because');
disp('of its size)');
disp(' ');
evalExacta=zeros(1,numPuntos);
for i=1:numPuntos
    evalExacta(i)=double(subs(r,sym(puntos(i))));
end

disp(' ');
disp('Computing the errors for VS and Casteljau algorithms, and Eval method');
ErrorAcc=abs((y_Eval(1,:)-evalExacta)./evalExacta);
ErrorVS=abs((y_VS-evalExacta)./evalExacta);
ErrorCast=abs((y_Cast-evalExacta)./evalExacta);
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp('Press any key to see a graphic with the errors');
pause

figure;
hold off;
semilogy(puntos,ErrorCast);
hold on;
semilogy(puntos,ErrorVS,'--');
semilogy(puntos,ErrorAcc,':');
hleg2 = legend('Casteljau','VS','Eval');
xlabel(texlabel('r(x)=x*(x-1/2)^20'));
ylabel('Relative errors');
xlim([0.4 0.6]);
hold off;

