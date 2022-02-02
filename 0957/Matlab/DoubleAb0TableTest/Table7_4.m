%TABLE7_4
%
global ab ab0
f0='%8.0f %5.1f %21.11e %12.4e\n';
f1='%14.1f %21.11e %12.4e\n';
dig=32; 
ab=loadvpa('ab_hrhermite',dig,200,2);
ab0=double(ab);

dig=20; digits(dig);
fprintf('\n')
disp('       n    x              y              err')        
errmax=0;
for n=[1:6 10 11]
  for x=0:.1:5
    if rem(x,.5)==0, fprintf('\n'), end

    y0  =  inerfc(n,x);                     % compute in double
    sy0 = sinerfc(dig,n,x);                 % compute in symbolic
    y   =     2^n *gamma(    1+n/2 )* y0;   % compute in double
    sy  = vpa(2^n)*gamma(vpa(1+n/2))*sy0;   % compute in symbolic

    err=double(abs((y-sy)/sy));             % compare double vs symbolic
    if err>errmax, errmax=err; end
    if x==0 
      fprintf(f0,n,x,y,err)
    else 
      fprintf(f1,x,y,err)
    end
  end
end
errmax
