%BOUND_FN_NEG
%
n=(0:150)';
y=zeros(151,1);
yl=log(2*(1+2.^n).*18.1.^(n+1)./(sqrt(pi)*gamma(n+1))+ ...
  sqrt(2)./gamma((n+2)/2));
[ylmax,N]=max(yl) 
ymax=exp(ylmax)
plot(n,yl);set(gca,'FontSize',14)
axis([0 150 -80 40])
xlabel('n'); ylabel('bound')
