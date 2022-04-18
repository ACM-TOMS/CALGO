clc
clear

matlab = ver('MATLAB');
matlabVersion = regexp(matlab.Version,'\.','split');
majorVersion = str2num(matlabVersion{1});
minorVersion = str2num(matlabVersion{2});
if majorVersion<7 || (majorVersion==7 && minorVersion<=1)
    error('Version of Matlab must be 7.2 at least');
end

n=8;

a1=0.24995;
b1=0.25005;
a2=0.74995;
b2=0.75005;
numPuntos=400;
puntos1=linspace(a1,b1,numPuntos);
puntos2=linspace(a2,b2,numPuntos);

x=sym('x','real');
p=sym((x-3/4)^7*(x-1));
q=sym(x*(x-1/4)^7);

disp(' ');
disp('This test file shows the performance of the Polynomial_1.0 package for the evaluation');
disp(' ');
disp('POLYNOMIALS considered: p(x)=(x-3/4)^7*(x-1) and q(x)=q(x)=x*(x-1/4)^7');
disp(' ')
disp(' ')
disp(' ');
disp('Computing the exact coefficients of p(x)=(x-3/4)^7*(x-1) with');
disp('respect to the Bernstein basis');

aux=sym(zeros(1,n+1));
for i=1:n+1
    aux(i)=sym(nchoosek(n,i-1))*x^(i-1)*(1-x)^(n+1-i);
end

c = sym('c', [1 n+1]);
if (majorVersion==7 && minorVersion >= 14) || majorVersion>=8
    exp1=taylor(expand(sum(c.*aux)),'Order',n+1);
    rhs=coeffs(taylor(expand(p),'Order',n+1),x);
else
    exp1=taylor(expand(sum(c.*aux)),n+1);
    rhs=coeffs(taylor(expand(p),n+1),x);
end

lhs=coeffs(exp1,x);
ecuacion=rhs-lhs;
ExactCoeffBer=solve(ecuacion);

for i=1:n+1
    coeffBer2(i)=double(eval(sprintf('ExactCoeffBer.c%d',i)));
    coeffVS2(i)=nchoosek(n,i-1)*coeffBer2(i);
end



disp(' ');
disp('Computing the exact coefficients of q(x)=x*(x-1/4)^7 with');
disp('respect to the Bernstein basis');

if (majorVersion==7 && minorVersion >= 14) || majorVersion>=8
    rhs=coeffs(taylor(expand(q),'Order',n+1),x);
else
    rhs=coeffs(taylor(expand(q),n+1),x);
end
rhs=[0,rhs];
ecuacion=rhs-lhs;
ExactCoeffBer=solve(ecuacion);

for i=1:n+1
    coeffBer1(i)=double(eval(sprintf('ExactCoeffBer.c%d',i)));
    coeffVS1(i)=nchoosek(n,i-1)*coeffBer1(i);
end

disp(' ');
disp('Evaluating q(x)=x*(x-1/4)^7 by VS algorithm at 400 points equally');
disp('distributed between 0.24995 and 0.25005, both included');
disp('(results are not shown in the screen because of its size)');

[y_VS1,errBound_VS1] = Vs(coeffVS1,puntos1);

disp(' ');
disp('Evaluating q(x)=x*(x-1/4)^7 by Casteljau algorithm at 400 points');
disp('equally distributed between 0.24995 and 0.25005, both included');
disp('(results are not shown in the screen because of its size)');

[y_Cast1,errBound_Cast1] = Casteljau(coeffBer1,puntos1);

disp(' ');
disp('Constructing the polynomial q(x)=x*(x-1/4)^7 from its exact coefficients respect'); 
disp('to the Bernstein basis with the corresponding constructor of the class');
disp('(Polynomial(coef,''c''))');

poly1 = Polynomial(coeffBer1,'c');

disp(' ');
disp('Evaluating q(x)=x*(x-1/4)^7 by Eval method at 400 points equally');
disp('distributed between 0.24995 and 0.25005 with e pretended precision');
disp('of the epsilon of the machine');
disp('(results are not shown in the screen because of its size)');

y_Eval1 = poly1.Eval(puntos1,1e-16);

disp(' ');
disp('Computing the exact values of the polynomial q(x) by using symbolic');
disp('capacities of Matlab (results are not shown in the screen because');
disp('of its size)');
disp(' ');
evalExacta1=zeros(1,numPuntos);
for i=1:numPuntos
    evalExacta1(i)=double(subs(q,sym(puntos1(i))));
end

disp(' ');
disp('Computing the errors for VS and Casteljau algorithms, and Eval method');
ErrorAcc1=abs(y_Eval1(1,:)-evalExacta1);
ErrorVS1=abs(y_VS1-evalExacta1);
ErrorCast1=abs(y_Cast1-evalExacta1);
ErrorAcc1_R=abs(y_Eval1(1,:)-evalExacta1)./abs(evalExacta1);
ErrorVS1_R=abs(y_VS1-evalExacta1)./abs(evalExacta1);
ErrorCast1_R=abs(y_Cast1-evalExacta1)./abs(evalExacta1);
disp(' ');

disp('Mean of the absolute errors when evaluating q(x) with the VS algorithm');
disp(mean(ErrorVS1));
disp('Mean of the absolute errors when evaluating q(x) with the de Casteljau algorithm');
disp(mean(ErrorCast1));
disp('Mean of the absolute errors when evaluating q(x) with the Eval method of the class');
disp(mean(ErrorAcc1));

disp(' ');
disp(' ');
disp(' ');
disp('Press any key to continue');
pause






disp(' ');
disp('Evaluating p(x)=(x-3/4)^7*(x-1) by VS algorithm at 400 points equally');
disp('distributed between 0.74995 and 0.75005, both included');
disp('(results are not shown in the screen because of its size)');

[y_VS2,errBound_VS2] = Vs(coeffVS2,puntos2);

disp(' ');
disp('Evaluating p(x)=(x-3/4)^7*(x-1) by Casteljau algorithm at 400 points');
disp('equally distributed between 0.74995 and 0.75005, both included');
disp('(results are not shown in the screen because of its size)');

[y_Cast2,errBound_Cast2] = Casteljau(coeffBer2,puntos2);

disp(' ');
disp('Constructing the polynomial p(x)=(x-3/4)^7*(x-1) from its exact coefficients respect'); 
disp('to the Bernstein basis with the corresponding constructor of the class');
disp('(Polynomial(coef,''c''))');

poly2 = Polynomial(coeffBer2,'c');

disp(' ');
disp('Evaluating p(x)=(x-3/4)^7*(x-1) by Eval method at 400 points equally');
disp('distributed between 0.74995 and 0.75005 with e pretended precision');
disp('of the epsilon of the machine');
disp('(results are not shown in the screen because of its size)');

y_Eval2 = poly2.Eval(puntos2,1e-16);

disp(' ');
disp('Computing the exact values of the polynomial p(x) by using symbolic');
disp('capacities of Matlab (results are not shown in the screen because');
disp('of its size)');
disp(' ');
evalExacta2=zeros(1,numPuntos);
for i=1:numPuntos
    evalExacta2(i)=double(subs(p,sym(puntos2(i))));
end

disp(' ');
disp('Computing the errors for VS and Casteljau algorithms, and Eval method');
ErrorAcc2=abs(y_Eval2(1,:)-evalExacta2);
ErrorVS2=abs(y_VS2-evalExacta2);
ErrorCast2=abs(y_Cast2-evalExacta2);
ErrorAcc2_R=abs(y_Eval2(1,:)-evalExacta2)./abs(evalExacta2);
ErrorVS2_R=abs(y_VS2-evalExacta2)./abs(evalExacta2);
ErrorCast2_R=abs(y_Cast2-evalExacta2)./abs(evalExacta2);
disp(' ');

disp('Mean of the absolute errors when evaluating q(x) with the VS algorithm');
disp(mean(ErrorVS2));
disp('Mean of the absolute errors when evaluating q(x) with the de Casteljau algorithm');
disp(mean(ErrorCast2));
disp('Mean of the absolute errors when evaluating q(x) with the Eval method of the class');
disp(mean(ErrorAcc2));

disp(' ')
disp('We can observe that Eval method is more accurate than VS and de Casteljau algorithms');
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp('Press any key to see a graphic with the errors');
pause

figure;
hold off;
subplot(1,2,1);
semilogy(puntos2,ErrorCast2_R);
hold on;
semilogy(puntos2,ErrorVS2_R,'--');
semilogy(puntos2,ErrorAcc2_R,':');
hleg2 = legend('Casteljau','VS','Eval');
xlabel(texlabel('p(x)=(x-3/4)^{7}(1-x)'));
ylabel('Relative errors');

subplot(1,2,2);
semilogy(puntos1,ErrorCast1_R);
hold on;
semilogy(puntos1,ErrorVS1_R,'--');
semilogy(puntos1,ErrorAcc1_R,':');
hleg1 = legend('Casteljau','VS','Eval');
xlabel(texlabel('p(t)=(x-1/4)^{7}x'));
ylabel('Relative errors');
hold off;
