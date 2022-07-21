function [innout,p_prime,alpha_coe,gap,iter_num]=Spherical_TA11(input_a,epsilon,alpha_0,query)
%%
%This function solves the convex hull membership query using Spherical
%Triangle algirithm. 
%
%Input:
%input_a: A mxn dimensional matrix representing set of points; m: dimension of set of points; n:
%numbe of points 
%epsilon: Scalar between (0,1); precision parameter.
%alpha_0: nx1 vector with non negative entries and sum up to 1; initialization of coefficient 
%query: mx1 vector; Query point in dimension m.
%Output:
%innout: {0,1} scalar. 0 if a witness is found, 1 if an epsilon apprixmation is found. 
%p_prime: mx1 vector. Witness of  epsilon approximation of query point.
%gap: euclidean distance between p_prime and query.
%iter_num: Number of iterations for algorithm to stop.
%%
[m,n]=size(input_a);
if nargin <=3
    mat_da=input_a;
    mat_norm_vec2=sum(mat_da.^2,1);
    matA_norm_vec=sqrt(mat_norm_vec2)+0.0001;
    veps=max(max(matA_norm_vec)*epsilon,epsilon);
    original_query = query;
    scale_matA=mat_da./matA_norm_vec;
else
    mat_da=input_a-repmat(query,1,n);
    mat_norm_vec2=sum(mat_da.^2,1);
    matA_norm_vec=sqrt(mat_norm_vec2)+0.000001;
    veps=max(max(matA_norm_vec)*epsilon,epsilon);
    original_query = query;
    scale_matA=mat_da./matA_norm_vec;
end
[m,n]=size(scale_matA);
artificial_data_n = round(0.01*n)+1;
rnd_A_raw=random('unif',0,1,n,artificial_data_n); 
norm_mat_A_raw=sum(rnd_A_raw,1);
rnd_A= rnd_A_raw./norm_mat_A_raw;
pts_A=scale_matA*rnd_A;
start_idx = n+1;
scale_matA = [scale_matA,pts_A];

norm_aug2=sum(scale_matA.^2,1);
[m,n]=size(scale_matA);
eudis=matA_norm_vec;
tmparr=find(eudis==min(eudis));
min_index=tmparr(1);
p_p=scale_matA(:,min_index);
alpha=zeros(1,n);
alpha(min_index)=1;
if(n<=1)
    innout=0;
    p_prime=p_p;
    alpha_coe=alpha;
    gap=min(eudis);
    return;
end

if(length(alpha_0)==n && n>1 && sum(alpha_0)==1)
    [ral,cal]=size(alpha_0);
    if ral>cal
        alpha_0=alpha_0';
    end
    alpha(1:n)=alpha_0;
    p_p=scale_matA*alpha';
end

innout=1;

beta_list=zeros(n,1);
iter_num=0;
distance=99999;
ppv=norm_aug2;
sp_dis=norm(p_p);
while(sp_dis>epsilon && distance>veps)
    iter_num=iter_num+1;
    gd=scale_matA'*(-p_p);
    norm2_ppv=(ppv+2*gd')+sum(p_p.^2);
    norm_ppv=sqrt(norm2_ppv);
    p_p_norm=sp_dis^2;
    dist_diff=-(p_p_norm+2*gd);
    index_pivot=find(dist_diff<0, 1);
    if isempty(index_pivot)
        found=0;
    else
        found=1;
        angle_val = p_p_norm+gd;
        pre_beta_list=angle_val./norm2_ppv';
        geq_index=find(pre_beta_list>=0);
        beta_list(geq_index)=min(1,pre_beta_list(geq_index));
        leq_index=find(pre_beta_list<0);
        beta_list(leq_index)=max(alpha(leq_index)'./(alpha(leq_index)-1-0.00000001)',pre_beta_list(leq_index));
        this_len=(abs(beta_list).*norm_ppv');
        v_index=find(this_len==max(this_len), 1 );
        beta=beta_list(v_index);
        alpha=(1-beta)*alpha;
        alpha(v_index)=alpha(v_index)+beta;
        p_p=(1-beta)*p_p+beta*scale_matA(:,v_index);
        rt_alpha = (alpha(start_idx:end)*rnd_A')+alpha(1:(start_idx-1));
        alpha_coe_rr=rt_alpha./matA_norm_vec;
        sum_rr=sum(alpha_coe_rr);
        sp_dis=norm(p_p);
        distance=sp_dis/sum_rr;

        
    end
    if(found==0)
        innout=0;
        break;
    end
end
rt_alpha=(alpha(start_idx:end)*rnd_A')+alpha(1:(start_idx-1));
alpha_coe_rr=rt_alpha./matA_norm_vec;
alpha_coe=alpha_coe_rr./(sum(alpha_coe_rr)+0.0000001);
p_prime=input_a*alpha_coe';
gap=norm(p_prime-original_query);
