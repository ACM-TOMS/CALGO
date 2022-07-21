function [p_prime,alpha_coe,dist,iter_num]=meb_algo(input_a,epsilon,query)
% 
[m,n]=size(input_a);
% disp('val epsilon1')
% epsilon
if nargin <=2
    mat_da=input_a;
    matA_norm_vec=sqrt(sum(mat_da.^2,1))+0.0001;
    scale_dmat= spdiags(matA_norm_vec(:).^(-1),0,n,n);
    scale_matA=mat_da*scale_dmat;
    veps=max(epsilon.*max(matA_norm_vec),epsilon);
    R_val=max(matA_norm_vec);
    scale_alpha=matA_norm_vec.^(-1);
    p=zeros(m,1);
    query=p;
else
    mat_da=input_a-query(:,ones(1,n));
    matA_norm_vec=sqrt(sum(mat_da.^2,1))+0.0001;
    scale_dmat= spdiags(matA_norm_vec(:).^(-1),0,n,n);
    scale_matA=mat_da*scale_dmat;
    veps=max(epsilon.*max(matA_norm_vec),epsilon);
    R_val=max(matA_norm_vec);
    scale_alpha=matA_norm_vec.^(-1);
    original_query = query;
    p=zeros(m,1);
    query=p;
    
end


[m,n]=size(scale_matA);



eudis=matA_norm_vec;

tmparr=find(eudis==min(eudis));
min_index=tmparr(1);
p_p=scale_matA(:,min_index);
alpha=zeros(1,n);
alpha(min_index)=1;
%disp('ok_3')
distance=min(eudis);
if(n<=1)
    inorout=0;
    p_prime=p_p;
    alpha_coe=alpha;
    dist=min(eudis);
    return;
end



pts_diff=p_p(:,ones(1,n))-scale_matA;
all_dis = sqrt(sum(pts_diff.^2,1));
max_dis = max(all_dis);
norm_v = sum(scale_matA.^2,1);
argmax_val = find(all_dis==max_dis,1);
delta = 0.5*max_dis;
r = 0.5*max_dis;
iter_num  = 0;
while delta> epsilon^2 && norm(p_p)>epsilon
%     epsilon
%     delta
    
    if iter_num>min(1/epsilon^2*3,100000)
        break;
    end
    for i=1:ceil(1/delta)
        iter_num=iter_num+1;
        beta = 1-r/max_dis;
        alpha = alpha*(1-beta);
        alpha(argmax_val) = alpha(argmax_val)+beta;
        p_p = (1-beta)*p_p + beta*scale_matA(:,argmax_val);
        inner_prod = -1*scale_matA'*p_p;
        pts_diff=  2*inner_prod'+ sum(p_p.^2) + norm_v;
        all_dis = sqrt(pts_diff);
        max_dis = max(all_dis);
        argmax_val = find(all_dis==max_dis,1);
        if norm(p_p)<epsilon
            break;
        end
    end
%     beta_val = (max_dis-1)
    s_val = max_dis -r;
    if s_val<=delta*3/4
        delta = delta*3/4;
    else
       r=r+delta/4;
       delta = delta*3/4;
    end
    
end
norm(p_p)
%     sqrt((p-p_p)'*(p-p_p))
% 
% dist=distance;
% % scale_alpha
% % alpha
alpha_coe_rr=alpha.*scale_alpha;
% % alpha_coe_rr
alpha_coe=alpha_coe_rr./(sum(alpha_coe_rr)+0.00000001);
dist=sqrt(sum((mat_da*alpha_coe'-query).^2));
p_prime=input_a*alpha_coe';
    
 