function []=test()
% M. Mehra and K. S. Patel %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%Input: It requires%%%%%%%%%%%%%%%
% dimension= 1 for one dimension and 2 for two dimension
% boundary condition type: boun_cond= 1 for Periodic and 2 for Dirichlet
% n= number of grid points
% x_l= left end in x direction
% x_r= right end in x direction
% y_l= left end in y direction
% y_r= right end in y direction
% f= function to be differentiated
% order(p)= Order of accracy desired
%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%
% Error= Error between the exact differentiation and compact finite
% difference approximation for a given function
dimension=input('Enter: 1 for one dimension, 2 for two dimesion: dimension=');
%%%%%%%%%%%%%%%%%   One Dimension %%%%%%%%%%%%%%%%%%%%%%%%%%
if dimension==1
 boun_type=input('Enter: 1 for Periodic Boundary, 2 for Dirichlet Boundary: boun_type=');
if boun_type==1
n=input('Enter the No. of grid points (Minimum 8): n=');
x_l=input('Enter the x min: x_l=');
x_r=input('Enter the x max: x_r=');
order=input('Enter: 1 for first derivative, 2 for second derivative, 3 for third derivative, 4 for fourth derivative: order=');
p=input('Enter the order of accuracy desired: 4 or 6 or 8 or 10: p=');
if order ==3 && p ~=4
    disp('If order is 3, Permisible value of p is only 6:')
    return
end
if order ==4 && p ~=4 && p ~=6
    disp('If order is 4, Permisible value of p are 4 and 6:')
    return
end
syms x
f=input('Enter a periodic function in variable x: f(x)=');
dx=(x_r-x_l)/n;y=x_l:dx:x_r-dx;fun=zeros(n,1);ana_df=zeros(n,1);
for i=1:n
    fun(i)=vpa(subs(f,x,y(i)));
end
if order==1
    [D]=compact_first_periodic(n,p,x_l,x_r);df=diff(f);size(D)
elseif order==2
    [D]=compact_second_periodic(n,p,x_l,x_r);df=diff(diff(f));
elseif order==3
    [D]=compact_third_periodic(n,p,x_l,x_r);df=diff(diff(diff(f)));
else
    [D]=compact_fourth_periodic(n,p,x_l,x_r);df=diff(diff(diff(diff(f))));
end
 deri=D*fun;
for i=1:n
ana_df(i)=vpa(subs(df,x,y(i)));
end
error=norm(ana_df-deri,inf);
disp('Error between analytic and compact finite difference approximation: Error=')
disp(error)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if boun_type==2
n=input('Enter the No. of grid points (Minimum 8): n=');
x_l=input('Enter the x min: x_l=');
x_r=input('Enter the x max: x_r=');
order=input('Enter: 1 for first derivative, 2 for second derivative: order=');
p=input('Enter the order of accuracy desired: 4 or 6: p=');
syms x
f=input('Enter the function in variable x: f(x)=');
y=linspace(x_l,x_r,n);fun=zeros(n,1);ana_df=zeros(n,1);
for i=1:n
    fun(i)=vpa(subs(f,x,y(i)));
end
if order==1
    [D]=first_compact_dirichlet(n,p,x_l,x_r);df=diff(f);
elseif order==2
    [D]=second_compact_dirichlet(n,p,x_l,x_r);df=diff(diff(f));
end
 deri=D*fun;
for i=1:n
ana_df(i)=vpa(subs(df,x,y(i)));
end
error=norm(ana_df-deri,inf);
disp('Error between analytic and compact difference approximation: Error=')
disp(error)
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%  Two Dimension  %%%%%%%%%%%%%%%
if dimension==2
 boun_type=input('Enter: 1 for Periodic boundary condition, 2 for Dirichlet boundary condition: boun_type=');
if boun_type==1
n=input('Enter the No. of grid points (Minimum 8): n=');
x_l=input('Enter the x min: x_l=');
x_r=input('Enter the x max: x_r=');
y_l=input('Enter the y min: y_l=');
y_r=input('Enter the y max: y_r=');
order=input('Enter 1 for d/dx, 2 for d/dy, 3 for d/dxx, 4 for d/dyy, 5 for d/dxy: order=');
p=input('Enter the order of accuracy desired: 4 or 6 or 8 or 10: p=');
syms x y
f=input('Enter the periodic function in variable x and y: f(x,y)=');
dx=(x_r-x_l)/n;dy=(y_r-y_l)/n;z=x_l:dx:x_r-dx;w=y_l:dy:y_r-dy;fun=zeros(n^2,1);ana_df=zeros(n^2,1);
a=1;
for j=1:n
for i=1:n
    fun(a)=vpa(subs(f,[x,y],[z(i),w(j)]));
    a=a+1;
end
end
if order==1
    [D]=compact_first_periodic_2dx(n,p,x_l,x_r);df=diff(f,x);
elseif order==2
    [D]=compact_first_periodic_2dy(n,p,y_l,y_r);df=diff(f,y);
elseif order==3
    [D]=compact_second_periodic_2dxx(n,p,x_l,x_r);df=diff(diff(f,x),x);
elseif order==4
    [D]=compact_second_periodic_2dyy(n,p,y_l,y_r);df=diff(diff(f,y),y);
else 
    [D]=compact_mixed_periodic_2dxy(n,p,x_l,x_r,y_l,y_r);df=diff(diff(f,x),y);
end
 deri=D*fun;
 b=1;
for j=1:n
    for i=1:n
ana_df(b)=vpa(subs(df,[x,y],[z(i),w(j)]));
b=b+1;
    end
end
error=norm(ana_df-deri,inf);
disp('Error between analytic and compact finite difference approximation: Error=')
disp(error)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if boun_type==2
n=input('Enter the No. of grid points (Minimum 8): n=');
x_l=input('Enter the x min: x_l=');
x_r=input('Enter the x max: x_r=');
y_l=input('Enter the y min: y_l=');
y_r=input('Enter the y max: y_r=');
order=input('Enter: 1 for d/dx, 2 for d/dy, 3 for d/dxx, 4 for d/dyy, 5 for d/dxy: order=');
p=input('Enter the order of accuracy desired: 4 or 6: p=');
syms x y
f=input('Enter the function in variable x and y: f(x,y)=');
z=linspace(x_l,x_r,n);w=linspace(y_l,y_r,n);fun=zeros(n^2,1);ana_df=zeros(n^2,1);
a=1;
for j=1:n
for i=1:n
    fun(a)=vpa(subs(f,[x,y],[z(i),w(j)]));
    a=a+1;
end
end
if order==1
    [D]=first_compact_dirichlet_2dx(n,p,x_l,x_r);df=diff(f,x);
elseif order==2
    [D]=first_compact_dirichlet_2dy(n,p,y_l,y_r);df=diff(f,y);
elseif order==3
    [D]=second_compact_dirichlet_2dxx(n,p,x_l,x_r);df=diff(diff(f,x),x);
elseif order==4
    [D]=second_compact_dirichlet_2dyy(n,p,y_l,y_r);df=diff(diff(f,y),y);
else 
    [D]=mixed_compact_dirichlet_2dxy(n,p,x_l,x_r,y_l,y_r);df=diff(diff(f,x),y);
end
 deri=D*fun;
 b=1;
for j=1:n
    for i=1:n
ana_df(b)=vpa(subs(df,[x,y],[z(i),w(j)]));
b=b+1;
    end
end
error=norm(ana_df-deri,inf);
disp('Error between analytic and compact finite difference approximation: Error=')
disp(error)
end
end
end