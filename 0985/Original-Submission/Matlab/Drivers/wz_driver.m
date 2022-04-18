clear all;close all;

ttiny = 0.06447*eps;
for icase=1:4
if icase==1
  y=logspace(-5,5,71);
  x=linspace(-500,500,40001);
elseif icase==2
  y=logspace(-20,4,71);
  x=linspace(-200,200,40001);  
elseif icase==3
    y=logspace(-5,5,71);
  x=linspace(-10,10,40001);    
elseif icase==4
clear x y
    y=logspace(-20,log10(6),71);
 for i=1:71     
     for j=1:40001
x(j)=-6+rand*2.*sqrt(36-y(i)^2);
 z(i,j)=x(j)+1i*y(i);
     end
 end
end
if icase ~=4
[x,y]=meshgrid(x,y);
z=complex(x,y);
end 

icase=icase

f_Hm=@()Humlicek_rev(z);         
V_Hm_time=timeit(f_Hm)/numel(z)


f_wz_lf=@()wz(z);         
V_wz_lf_time=timeit(f_wz_lf)/numel(z)

Ratio=V_wz_lf_time/V_Hm_time
end