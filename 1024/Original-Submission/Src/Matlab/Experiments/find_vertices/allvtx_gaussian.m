Num_rep=10;
Dim_array=[5,7,9,10,100,200];
Vtx_array=[50,100,200,300,500,1000];

frac_array=[0,20/100,50/100];
time_ta=zeros(length(frac_array),length(Dim_array));
time_scale_ta=zeros(length(frac_array),length(Dim_array));
time_qh=zeros(length(frac_array),length(Dim_array));
problem_size=zeros(1,length(frac_array));
for ii=1:length(Dim_array)
    Num_dim=Dim_array(ii);
    

    for jj=1:length(frac_array)
        frac_rat=frac_array(jj);
        Num_vtx=ceil(Vtx_array(ii)*(1-frac_rat))
        Num_pts=ceil(Vtx_array(ii)*(frac_rat))
        for kk=1:Num_rep
            ii
            jj
            kk
            vtx_A=Random_pts(Num_dim,Num_vtx,'normal');
            if frac_rat~=0
                inhulldata=Random_cvx(vtx_A,Num_pts,'dir');
                matA=[inhulldata,vtx_A];
                [m,n]=size(matA);
                normal_val_index=(Num_pts+1):n;

                rndindx=randperm(n,n);
                [or_val or_ind]=sort(rndindx);
                vertices_index=or_ind(normal_val_index);
                matA=matA(:,rndindx);
            else
                matA=vtx_A;
                [m,n]=size(matA);
                
            end
            disp('ta')
            tic;
            [index_this_avt,iter_ta_val]=AVTA_anti(matA',0.0001);
            ta_end=toc;
            time_ta(jj,ii)=time_ta(jj,ii)+ta_end;
            disp('end ta')
            disp('spta')
            tic;
            [index_this_spavt,iter_spta]=AVTA_SP(matA',0.0001);
            scale_ta_end=toc;
            time_scale_ta(jj,ii)=time_scale_ta(jj,ii)+scale_ta_end;
            disp('end spta')
            

            
            
            
            disp('qh')
            if ii<=4
                [QH,v] = convhulln(matA');
                qh_t=toc;
                time_qh(jj,ii)=time_qh(jj,ii)+qh_t;
                disp('end qh')
            end
%             
            
        end
        
    end
    
    
end
        




names=strings(6,1);
for i=1:length((time_ta(1,:)))
    this_size=[num2str(Dim_array(i)), 'x', num2str(Vtx_array(i))];
    names(i)=this_size;
end
        
for i=1:length((time_ta(:,1)))

    this_title=['vertice generated using Gaussian', num2str(frac_array(i)*100),'% redundant point'];
    figure(i);
    hold on
%     title(this_title);

    plot(log(time_ta(i,:)/Num_rep),'DisplayName','AVTA','LineWidth',1.5)
    plot(log(time_scale_ta(i,:)/Num_rep),'DisplayName','AVTA+','LineWidth',1.5)  
    plot(log(time_qh(i,:)/Num_rep),'DisplayName','Quick Hull','LineWidth',1.5)  ;
    legend('show','Location','northwest','fontsize',10)%,'Orientation','horizontal')
    set(gca,'xtick',[1:length((time_ta(i,:)))],'xticklabel',names,'fontsize',15)
    xlabel ("Problem Size",'fontsize',15);
    ylabel ("Running time (secs in log scale)",'FontSize',15);
    save_title= strrep(this_title,' ','_');
    saveas(gcf,['allvertices_',save_title,'.png'])
    hold off;
end

