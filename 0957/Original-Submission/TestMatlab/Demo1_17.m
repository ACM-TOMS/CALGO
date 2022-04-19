%DEMO1.17
%
f0='%2.0f %22.14e %22.14e %11.4e %11.4e\n';
disp(' k            x                     w                errx       errw')
dig=32;
sab=loadvpa('ab_hrhermite',3,dig,200,2);
ab=double(sab); 
N=40;
xw=gauss(N,ab); sxw=sgauss(dig,N,sab);
%sw=sxw(:,2)
errx=double(abs((xw(:,1)-sxw(:,1))./sxw(:,1)));
errw=double(abs((xw(:,2)-sxw(:,2))./sxw(:,2)));
for k=1:N
  fprintf(f0,k,xw(k,1),xw(k,2),errx(k),errw(k))
end
