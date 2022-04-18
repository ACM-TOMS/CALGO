clear all
clc

%% Point of Evaluation %%

a=0.5;

%% Function %%

syms x

% g=exp(x)*((x)^2);
% g=(exp(x)^2)+x+(exp(x));
% g=(x+exp(x))^0.5;
 g=((sin(x)+((x^2)/cos(x)))^0.5);
% g=(((x)^2)+exp(x))^0.5;
% g=sqrt((exp(x)^2)+x+(exp(x)));
% g=log(acos(x))+x;

y1=diff(g,1);
y2=diff(g,2);
y3=diff(g,3);
y4=diff(g,4);
y5=diff(g,5);
y6=diff(g,6);
y7=diff(g,7);


D1=vpa(subs(y1,x,a));
D2=vpa(subs(y2,x,a));
D3=vpa(subs(y3,x,a));
D4=vpa(subs(y4,x,a));
D5=vpa(subs(y5,x,a));
D6=vpa(subs(y6,x,a));
D7=vpa(subs(y7,x,a));


% D2=eval(subs(y,x,a));


% G = @(x) exp(x)*((x)^2);
% G = @(x) (exp(x)^2)+x+(exp(x));
% G = @(x) (x+exp(x))^0.5;
G = @(x) ((sin(x)+((x^2)/cos(x)))^0.5);
% G = @(x) (((x)^2)+exp(x))^0.5;
% G = @(x) ((exp(x)^2)+x+(exp(x)))^0.5;
% G = @(x) log(acos(x))+x;


for j=1:15
    
h(j)=10^(-j);


%% Complex Step Multicomplex Toolbox %%

D1F=(((G(multicomplex(inputconverter(a,1,h(j))))))/(h(j)));
errormulti1(j)=abs((D1F.zn(end)-D1)*100/D1);
D2F=(((G(multicomplex(inputconverter(a,[1,2],h(j))))))/(h(j)^2));
errormulti2(j)=abs((D2F.zn(end)-D2)*100/D2);
D3F=(((G(multicomplex(inputconverter(a,[1,2,3],h(j))))))/(h(j)^3));
errormulti3(j)=abs((D3F.zn(end)-D3)*100/D3);
D4F=(((G(multicomplex(inputconverter(a,[1,2,3,4],h(j))))))/(h(j)^4));
errormulti4(j)=abs((D4F.zn(end)-D4)*100/D4);
D5F=(((G(multicomplex(inputconverter(a,[1,2,3,4,5],h(j))))))/(h(j)^5));
errormulti5(j)=abs((D5F.zn(end)-D5)*100/D5);
D6F=(((G(multicomplex(inputconverter(a,[1,2,3,4,5,6],h(j))))))/(h(j)^6));
errormulti6(j)=abs((D6F.zn(end)-D6)*100/D6);
D7F=(((G(multicomplex(inputconverter(a,[1,2,3,4,5,6,7],h(j))))))/(h(j)^7));
errormulti7(j)=abs((D7F.zn(end)-D7)*100/D7);

disp('loop loop loop loop:')
j
end

figure(1) 

loglog(h,errormulti1,'+-')
hold on
loglog(h,errormulti2,'o-')
hold on
loglog(h,errormulti3,'^-')
hold on
loglog(h,errormulti4,'d-')
hold on
loglog(h,errormulti5,'x-')
hold on
loglog(h,errormulti6,'s-')
hold on
loglog(h,errormulti7,'v-')
% holdThird Derivative error for different h, $$g=\sqrt{e^{2x}+x+e^{x}}$$','FontSize',18,'Interpreter','Latex')
xlabel('$$h$$','FontSize',15,'Interpreter','Latex')
ylabel('$$error \enskip (\%) $$','FontSize',15,'Interpreter','Latex')
legend('1st Derivative','2nd Derivative','3rd Derivative','4th Derivative','5th Derivative','6th Derivative','7th Derivative','Location','northeast')
grid on
grid minor
% xlim([10^-15 10^0])
% ylim([11^-15 10^5])

set(gcf,'PaperPosition',[0 0 14 10.5])
print -dpng ErrorComparisonFigure

