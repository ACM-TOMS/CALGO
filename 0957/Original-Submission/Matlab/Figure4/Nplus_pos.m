%NPLUS_POS  The number of Gauss points needed when x>=0
%
global ab ab0
digits(32); dig=32;
ab=loadvpa('ab_hrhermite',dig,200,2);
ab0=double(ab); 
N0=5; dN=1; 
eps0=.5d-12;
n0=0;
%n0=25;
n0=50;
%n0=75;
%n0=100;
%n0=125;
i=0;
for x=0:.2:27
  if x<=14.4
    i=i+1; xp(i)=x;
    Nm=0;
    for n=n0:n0+25
      N=N_inerfc(n,x,eps0,N0,dN);
      if N>Nm, Nm=N; end
    end
    Np(i)=Nm;
  else
    nmax=floor((27-x)/.084);
    if nmax>=n0
      Nm=0; i=i+1; xp(i)=x;
      for n=n0:nmax
        [N,y]=N_inerfc(n,x,eps0,N0,dN);
        if N>Nm, Nm=N; end
      end
      Np(i)=Nm;
    end
  end
end  
%  [xp' Np']
plot(xp,Np);set(gca,'FontSize',14)
xlabel('x'); ylabel('N')
%  axis([0 20 34 56]) %only for n0=100
