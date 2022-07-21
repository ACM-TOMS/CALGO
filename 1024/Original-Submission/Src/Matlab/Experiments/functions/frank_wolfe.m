function [alpha err k]=frank_wolfe(p,A,epsilon)
 
 
%%*********Frank wolfe algorithm minimize eucidean distance(Update rule )*********
%Input:
%   p: mx1 vector The point which needs to be approximated
%   A: mxn matrix. Columns of A are points in convex hull.
%   epsilon: precision parameter
%Output:
%   alpha: nx1 vector. The weight of vertices to approximate p.
%   err: The l1-distance of p to the convex hull. 
%%***************************************************************
 
 
[m,n]=size(A);
 
alpha=zeros(n,1);
diff=A-p(:,ones(n,1));
dis=sum(diff.^2,1);
min_ind=min(find(dis==min(dis)));
max_dis = max(sqrt(dis));
alpha(min_ind)=1;
err=min(dis);
if(n<=1)
    return;
end
 
 
p_p=A(:,min_ind);
gd=A'*(p_p-p);
mu=-1*min(gd);
lambda=gd+mu;
 
 
 
k=0;
 
while k<min(100/epsilon,1000000)
    k=k+1;
    if(sqrt(sum((p-p_p).^2))<epsilon*max_dis)     
       break; 
    end
    s=min(find(gd==min(gd)));
    gamma=2/(k+3);
    alpha=(1-gamma).*alpha;
    alpha(s)=alpha(s)+gamma;
    p_p=(1-gamma)*p_p+gamma*A(:,s);
    gd1=A'*p_p;
    gd2=A'*p;
    gd=gd1-gd2;
end
 
err=sqrt(sum((p-p_p).^2));
 
    
