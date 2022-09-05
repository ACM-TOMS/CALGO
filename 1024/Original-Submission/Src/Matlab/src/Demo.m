%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a demo for Spherical Triangle Algrithm
%% Full paper:
%% Spherical Triangle Algorithm: A Fast 
%% Oracle for Convex Hull Membership Queries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clear all


Num_rep=10; %% Number of repeats
Problem_size=[100,1000;200,2000;500,5000;1000,10000]; %% Problem size (dim,number of points)
epsilon=0.001; %% Precision parameter

time_ta=zeros(1, length(Problem_size(:,1))); %%Triangle Algorithm running time
time_sphere_ta=zeros(1, length(Problem_size(:,1))); %%Spherical Triangle Algorithm running time
time_lp=zeros(1, length(Problem_size(:,1)));  %% Linear Programming running time
iter_ta=zeros(1, length(Problem_size(:,1)));     %% Itertion number Spherical Triangle Algorithm
iter_sphere_ta=zeros(1, length(Problem_size(:,1))); %% Itertion number Triangle Algorithm

for ii=1:length(Problem_size(:,1))
    Num_dim=Problem_size(ii,1);
    Num_vtx=Problem_size(ii,2);
    for kk=1:Num_rep
        %%%Data generation
        vtx_A=Random_pts(Num_dim,Num_vtx,'normal');
        matA=vtx_A;
        [m,n]=size(matA);
        p=Random_cvx(vtx_A,1,'dir');%% Query points
%         p=Random_pts(Num_dim,1,'normal');
        %%%
        disp('Triangle Algorithm starts')
	  
        %%Triangle Algorithm
        tic;
        [inorout,p_prime,alpha_coe,dist,ta_iter]=Tri_Algo(matA,epsilon,p);
        ta_end=toc;
        time_ta(ii)=time_ta(ii)+ta_end;
        iter_ta(ii)=iter_ta(ii)+ta_iter;
        disp('Triangle Algorithm time')
        disp(ta_end)
        %%Spherical Triangle Algorithm
        disp('Spherical Triangle Algorithm  starts')
        tic;
        [inorout,p_prime,alpha_coe,dist,iter_num]=Spherical_TA(matA,epsilon,zeros(n,1),p);
        sphere_ta_end=toc;
        time_sphere_ta(ii)=time_sphere_ta(ii)+sphere_ta_end;
        iter_sphere_ta(ii)=iter_sphere_ta(ii)+iter_num;
        disp('Spherical Triangle Algorithm  time')
        disp(sphere_ta_end)
        
        
        
        
        %%Linear Programming
        disp('Linear Programming starts')
        tic;

        n_var=n;
        n_con=m;
        x_y_size=n_con+n_var;
        cvx_aeq=ones(1,n_var);
        Aeq=[matA,eye(n_con);ones(1,n_var),zeros(1,n_con)];
        beq=[p;1];
        lb=zeros(x_y_size,1);
        tic
        options = optimset('linprog');
        options.Display = 'off';

        [bb,fval] = linprog([zeros(n_var,1);ones(n_con,1)]...
        ,[],[],Aeq,beq,lb,[],options);
        lp_t=toc;
        time_lp(ii)=time_lp(ii)+lp_t;
        disp('Linear Programming time')
        disp(lp_t)
    end  
end
        

hold on
% title('query point infeasible')
title('query point feasible')
plot(log((time_ta(1,:))),'DisplayName','TA','LineWidth',1.5)
plot(log((time_sphere_ta(1,:))),'DisplayName','Spherical TA','LineWidth',1.5)  
plot(log(time_lp(1,:)),'DisplayName','LP','LineWidth',1.5)  ;
legend('show','Location','northwest')%,'Orientation','horizontal')
hold off;
