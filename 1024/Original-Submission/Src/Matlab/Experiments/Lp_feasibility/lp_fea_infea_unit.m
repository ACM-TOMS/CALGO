Num_rep=10;
Dim_array=[50,100,200];
Vtx_array=[500,1000,2000];
Eps_array=[0.01,0.01,0.01];
time_ta=zeros(length(Eps_array),length(Dim_array));
time_scale_ta=zeros(length(Eps_array),length(Dim_array));
time_lp=zeros(length(Eps_array),length(Dim_array));
time_mebalgo=zeros(length(Eps_array),length(Dim_array));
time_fw=zeros(length(Eps_array),length(Dim_array));
time_eg=zeros(length(Eps_array),length(Dim_array));
iterrr=0;
for ii=1:length(Eps_array)
    Num_dim=Dim_array(ii);
    Num_vtx=Vtx_array(ii);
    
    for jj=1:length(Eps_array)
        epsilon=0.00001;
        for kk=1:Num_rep
            ii
            jj
            kk
            iterrr=iterrr+1;
            A=Random_pts(Num_dim,Num_vtx,'unit ball');
            [AL,Sig,AR] = svd(A);
            Sig(1:round(Num_dim/2),:)=0;
            X=Random_pts(Num_vtx,1,'normal');
            b = AL*(Sig+Random_pts(size(Sig,1),size(Sig,2),'normal'))*AR*X;
            A=AL*Sig*AR;
            [M_con,N_var]=size(A);
            
            M=1200;
            tmp_mat=([A,zeros(M_con,1);ones(1,N_var),1]);
            tmp_b=[-b;-M];
            data_mat=[tmp_mat,tmp_b;zeros(1,N_var+1),1];
            p=[zeros(Num_dim+1,1);1/(1+M)];
            disp('ta')
            tic;
            [inorout,p_prime_1,alpha_coe,dist,ta_iter]=ta_anti(data_mat,p,epsilon);
            ta_end=toc;
            time_ta(ii,jj)=time_ta(ii,jj)+ta_end;
            disp('end ta')
            disp('spta')
            tic;
            [inorout,p_prime_2,alpha_coe,dist,iter_num]=Spherical_TA11(data_mat,epsilon,[0,0],p);
            scale_ta_end=toc;
            time_scale_ta(ii,jj)=time_scale_ta(ii,jj)+scale_ta_end;
            
            disp('end spta')
            sta_x=alpha_coe(1:(end-2))/alpha_coe(end-1);
            
            
            
            disp('fw')
            tic;
            [alpha err iter_num]=frank_wolfe(p,data_mat,epsilon);
            fw_end=toc;
            time_fw(ii,jj)=time_fw(ii,jj)+fw_end;
            disp('end fw')
            
            disp('meb')
            tic;
            [p_p,alpha,argmax_val,iter_num]=meb_algo(data_mat,epsilon,p);
            mebalgo_end=toc;
            time_mebalgo(ii,jj)=time_mebalgo(ii,jj)+mebalgo_end;
            
            disp('end meb')
            
            
            disp('lp')
            tic;
            [n_1,n_2]=size(A);
            x_y_size=n_2;
            c_coe=[zeros(n_2,1)];

            Aeq=[A];

            beq=b;
            lb=zeros(x_y_size,1);
            
            [bb,fval] = linprog([],[A;-A],[b;-b],[],[],lb,[]);
            lp_t=toc;
            time_lp(ii,jj)=time_lp(ii,jj)+lp_t;
            disp('end lp')
            
            
            
        end
        
    end
    
    
end
        
        

names=strings(3,1);
for i=1:length((time_ta(:,1)))
    this_size=[num2str(Dim_array(i)), 'x', num2str(Vtx_array(i))];
    names(i)=this_size;
end
        
for i=1:length((time_ta(:,1)))
    this_title=['Unit sphere infeasible case'];
    figure(i);
     hold on
%     title(this_title);
    plot(log(time_ta(:,i)/Num_rep),'DisplayName','Triangle Algorithm','LineWidth',1.5)
    plot(log(time_scale_ta(:,i)/Num_rep),'DisplayName','Spherical TA','LineWidth',1.5)  
    plot(log(time_lp(:,i)/Num_rep),'DisplayName','Lp Solver','LineWidth',1.5)  ;
    plot(log(time_fw(:,i)/Num_rep),'DisplayName','Frank Wolfe','LineWidth',1.5)  ;
    plot(log(time_mebalgo(:,i)/Num_rep),'DisplayName','MEBOPT','LineWidth',1.5)  ;
    legend('show','Location','northwest','fontsize',10)%,'Orientation','horizontal')
    set(gca,'xtick',[1:length((time_ta(:,1)))],'xticklabel',names,'fontsize',15)
    xlabel ("Problem Size",'FontSize',20);
    ylabel ("Running time (secs in log scale)",'FontSize',15);
    save_title= strrep(this_title,' ','_');
    saveas(gcf,['lp_feasibility_',save_title,'.png'])
    hold off;
end


