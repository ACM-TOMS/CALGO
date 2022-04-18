%RUN_XUNDER
%
global ab ab0
digits(32); dig=32;
ab=loadvpa('ab_hrhermite',dig,200,2);
ab0=double(ab); 
x=zeros(151,1);
i=0;
for n=0:150
  i=i+1; x(i)=xunder(n);
end
n=(0:150)';
xa=27-12.6*n/150;
plot(n,x);set(gca,'FontSize',14)
xlabel('n'); ylabel('x')
hold on
plot(n,xa,'--')
