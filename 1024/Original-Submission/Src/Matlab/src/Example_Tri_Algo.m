m=10;
n=20;
epsilon=0.01;
matA=Random_pts(m,n,'unit ball'); %% Generate points from unit sphere
p=Random_cvx(matA,1);%% query points p is in the convex hull
tic
[inorout,p_prime,alpha_coe,gap,iter_num]...
=Tri_Algo(matA,epsilon,p);
toc
alpha_coe
p_prime 
inorout
gap