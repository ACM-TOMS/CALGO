function [inorout,p_prime,alpha_coe,dist,iter]=ta_anti(dataset,p,epsilon)



matA=dataset;

[m,n]=size(matA);





diffmat=matA-p(:,ones(1,n));
eudis=sqrt(sum(diffmat.^2,1));
veps=max(epsilon*max(eudis),epsilon);
tmparr=find(eudis==min(eudis));
min_index=tmparr(1);
p_p=matA(:,min_index);
alpha=zeros(1,n);
alpha(min_index)=1;
distance=min(eudis);
if(n<=1)
    inorout=0;
    p_prime=p_p;
    alpha_coe=alpha;
    dist=min(eudis);
    return;
end

iter=0;
inorout=1;
dist_vp=sum(diffmat.*diffmat,1);
beta_list=zeros(n,1);
while(sqrt((p-p_p)'*(p-p_p))>veps)
    iter=iter+1;
    found=0;
    distance=sqrt((p-p_p)'*(p-p_p));
    gd1=matA'*p;
    gd2=matA'*p_p;
    gd=gd1-gd2;
    ppv=(dist_vp+2*gd')+sum(p_p.^2)-p'*p;
    norm2_ppv=ppv;
    norm_ppv=sqrt(norm2_ppv);
    p_norm=p'*p;
    p_p_norm=p_p'*p_p;
    dist_diff=(p_norm-p_p_norm)- 2*gd;
    index_pivot=find(dist_diff<=0);
    pre_beta_list=(p_p'*(p_p-p)+gd)./norm2_ppv';
    if length(index_pivot)==0
        found=0;
    else
        geq_index=find(pre_beta_list>=0);
        beta_list(geq_index)=min(1,pre_beta_list(geq_index));
        leq_index=find(pre_beta_list<0);
        beta_list(leq_index)=max(alpha(leq_index)'./(alpha(leq_index)-1-0.0000001)',pre_beta_list(leq_index));
        this_len=abs(beta_list).*norm_ppv';
        v_index=min(find(this_len==max(this_len)));
        beta=beta_list(v_index);
        alpha=(1-beta)*alpha;
        alpha(v_index)=alpha(v_index)+beta;
        p_p=matA*alpha';
        found=1;
    end
    if(found==0)
        inorout=0;
        break;
    end
end

p_prime=p_p;
alpha_coe=alpha;
dist=distance;