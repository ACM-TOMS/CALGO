Num_rep=10;
Dim_array=[100,500,1000];
Vtx_array=[500,2000,5000];
Eps_array=[0.001,0.0001,0.00001];
time_ta=zeros(length(Eps_array),length(Dim_array));
time_scale_ta=zeros(length(Eps_array),length(Dim_array));
time_lp=zeros(length(Eps_array),length(Dim_array));
time_mebalgo=zeros(length(Eps_array),length(Dim_array));
time_fw=zeros(length(Eps_array),length(Dim_array));
time_eg=zeros(length(Eps_array),length(Dim_array));
for ii=1:length(Dim_array)
    Num_dim=Dim_array(ii);
    Num_vtx=Vtx_array(ii);

    for jj=1:length(Eps_array)
        epsilon=Eps_array(jj);
        for kk=1:Num_rep
            vtx_A=Random_pts(Num_dim,Num_vtx,'unit ball');
            matA=vtx_A;
            [m,n]=size(matA);
            p=Random_cvx(vtx_A,1,'dir');
            disp('ta')
            tic;
            [inorout,p_prime,alpha_coe,dist,ta_iter]=ta_anti(matA,p,epsilon);
            ta_end=toc;
            time_ta(ii,jj)=time_ta(ii,jj)+ta_end;
            disp('end ta')
            disp('spta')
            tic;
            [inorout,p_prime,alpha_coe,dist,iter_num]=Spherical_TA11(matA,epsilon,[0,0],p);
            scale_ta_end=toc;
            time_scale_ta(ii,jj)=time_scale_ta(ii,jj)+scale_ta_end;
            disp('end spta')
            tic;
            disp('meb')
            [p_p,alpha,argmax_val,iter_num]=meb_algo(matA,epsilon,p);
            mebalgo_end=toc;
            time_mebalgo(ii,jj)=time_mebalgo(ii,jj)+mebalgo_end;
            disp('end meb') 
            tic;
            disp('fw')
            [alpha err iter_num]=frank_wolfe(p,matA,epsilon);
            fw_end=toc;
            time_fw(ii,jj)=time_fw(ii,jj)+fw_end;
            disp('end fw')
            
            disp('lp')
            tic;
            
            n_var=n;
            n_con=m;
            x_y_size=n_con+n_var;
            cvx_aeq=ones(1,n_var);
            Aeq=[matA,eye(n_con);ones(1,n_var),zeros(1,n_con)];
            beq=[p;1];
            lb=zeros(x_y_size,1);
            tic
            [bb,fval] = linprog([zeros(n_var,1);ones(n_con,1)],[],[],Aeq,beq,lb,[]);
            lp_t=toc;
            time_lp(ii,jj)=time_lp(ii,jj)+lp_t;
            disp('end lp')
            
            
            
        end
        
    end
    
    
end
        

names=strings(4,1);
for i=1:length((time_ta(:,1)))
    this_size=[num2str(Dim_array(i)), 'x', num2str(Vtx_array(i))];
    names(i)=this_size;
end
        
for i=1:length((time_ta(:,1)))
    this_title=['Unit Sphere feasible case ','epsilon ',num2str(Eps_array(i))];
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
    saveas(gcf,[save_title,'.png'])
    hold off;
end










