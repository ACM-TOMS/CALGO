Num_rep=10;
Dim_array=[20,50,100];
Num_vtx=500;
Red_array=[2000,10000,20000];
Eps_array=[0.001,0.001,0.001];
time_ta=zeros(length(Eps_array),length(Dim_array));
time_scale_ta=zeros(length(Eps_array),length(Dim_array));
time_qh=zeros(length(Eps_array),length(Dim_array));


for ii=1:length(Dim_array)

    Num_dim=Dim_array(ii);
    
    epsilon=Eps_array(ii);
    for jj=1:length(Eps_array)
        Num_pts=Red_array(jj);
        
        for kk=1:Num_rep
            ii
            jj
            kk
            vtx_A=Random_pts(Num_dim,Num_vtx,'normal');
            inhulldata=Random_cvx(vtx_A,Num_pts,'dir');
            matA=[inhulldata,vtx_A];
            [m,n]=size(matA);
            normal_val_index=(Num_pts+1):n;

            rndindx=randperm(n,n);
            [or_val or_ind]=sort(rndindx);
            vertices_index=or_ind(normal_val_index);
            matA=matA(:,rndindx);
            
            disp('ta')
            tic;
            [index,iter_num_ta]=AVTA_anti(matA',0.0001);
            [avta_A , avta_c] = MinVolEllipse(matA(:,index),epsilon);
            disp('the running time of MVE with  AVTA is')
            ta_end=toc
            time_ta(jj,ii)=time_ta(jj,ii)+ta_end;
            disp('end ta')
            tic;
            [index iter_sp_ta]=AVTA_SP(matA',0.0001);
            [spavta_A , spavta_c] = MinVolEllipse(matA(:,index),epsilon);
            disp('the running time of MVE with spherical AVTA is')
            scale_ta_end=toc
            time_scale_ta(jj,ii)=time_scale_ta(jj,ii)+scale_ta_end;
            disp('end spta')
            disp('mvee')
            tic;
            [A , c] = MinVolEllipse(matA,epsilon);
            disp('the running time of MVE without AVTA is')
            qh_t=toc
            time_qh(jj,ii)=time_qh(jj,ii)+qh_t;
            disp('end mvee')
            
            
            
        end
        
    end
    
    
end

ylabels = string(Red_array);
xlabels = string(Dim_array);
x_name = 'Dimension';
y_name = 'Number of Redundant Points';


title_name = 'MVEE AVTA+ running time (secs)'
h=heatmap(time_scale_ta/Num_rep,'FontSize',17)


h.XDisplayLabels=xlabels
h.YDisplayLabels=ylabels
h.XLabel = 'Dimension'
h.YLabel = 'Number of Redundant Points'
% h.Title = title_name
h.caxis([0, 500]);
save_title= strrep(title_name,' ','_');
h.ColorbarVisible = 'off';
saveas(gcf,[save_title,'.png'])

title_name = 'MVEE AVTA running time (secs)'
h=heatmap(time_ta/Num_rep,'FontSize',17)
h.XDisplayLabels=xlabels
h.YDisplayLabels=ylabels
h.XLabel = 'Dimension'
% h.YLabel = 'Number of Redundant Points'
% h.Title = title_name
h.caxis([0, 500]);
h.ColorbarVisible = 'off';
save_title= strrep(title_name,' ','_');
saveas(gcf,[save_title,'.png'])


title_name = 'MVEE running time (secs)'
h=heatmap(time_qh/Num_rep,'FontSize',17)
h.XDisplayLabels=xlabels
h.YDisplayLabels=ylabels
h.XLabel = 'Dimension'
% h.YLabel = 'Number of Redundant Points'
% h.Title = title_name
h.caxis([0, 500]);
save_title= strrep(title_name,' ','_');
saveas(gcf,[save_title,'.png'])

        


