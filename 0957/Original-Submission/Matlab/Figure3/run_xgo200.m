%RUN_XGO200
%
global ab ab0
digits(32); dig=32;
ab=loadvpa('ab_hrhermite',dig,200,2);
ab0=double(ab); 
x=zeros(151,1);
i=0;
for n=0:150
  i=i+1; x(i)=xgo200(n);
end
n=(0:150)';
xa=-17.9+.024*n;
plot(n,x);set(gca,'FontSize',14)
xlabel('n'); ylabel('x')
hold on
plot(n,xa,'--')
