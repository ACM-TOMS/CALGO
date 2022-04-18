%NPLUS_NEG The number of Gauss point needed when x<0
% 
global ab ab0
f0='%8.0f %4.0f %7.4f\n';
digits(32); dig=32;
ab=loadvpa('ab_hrhermite',dig,200,2);
ab0=double(ab); 
N0=5; dN=1;;
eps0=.5e-12; 
nmax=150;
for x=-.9:-1:-16.9
  if x<-14.3, nmax=floor((x+17.9)/.024); end
  N=zeros(nmax,1);
  for n=1:nmax
    N(n)=N_inerfc(n,x,eps0,N0,dN);
  end
  N0=N(1); 
  [N0 (N(nmax)-N0)/(nmax-1)]
  n=(1:nmax)';
  plot(n,N);set(gca,'FontSize',14)
  axis('square')
  axis([0 210 0 210])
  xlabel('n'); ylabel('N')
  hold on
  text(160,60,'x = - 0.9','FontSize',14)
  text(160,90,'x = - 4.9','FontSize',14)
  text(160,130,'x = -8.9','FontSize',14)
  text(160,180,'x = -12.9','FontSize',14)
  text(25,204,'x = -16.9','FontSize',14)
end
