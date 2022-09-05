function [index]=AVTA_PLUS(matrix_A,GAMMA,varargin)

%%*********Compute K vertices using avta_K
%Input: 
%   matrix_A: dimension x sample number matrix. Rows are points.
%   GAMMA: The precision parameter
%   varargin: 'Rp', dim of projection,  if using random projection, 
%Output:
%   index: index of vertices
%%***************************************
matrix_A=matrix_A';
Type='';
if ~isempty(varargin)
    if nargin>=4
        [Type Para]=deal(varargin{:});
    else
        Type=deal(varargin{:});
        Para=size(matrix_A,2)*2;
    end
end

if strcmp(Type,'Rp')==1
    [NN MM]=size(matrix_A);
    rand_matrix=random('Normal',0,1,MM,Para);
    matrix_A_raw=matrix_A;
    matrix_A=matrix_A*rand_matrix;
end



[NN,MM]=size(matrix_A);

index_set=[];


epsilon=2;

eudis=sum(matrix_A.*matrix_A,2);
tmparr=find(eudis==max(eudis));
max_index=tmparr(1);


index_set=logical(zeros(1,NN));
vtx_ind_1=max_index;

index_set(vtx_ind_1)=1;


S_set=logical(ones(1,NN));
candidate_set=S_set;
candidate_set(index_set)=0;
index_series=1:NN;
alpha_zero=ones(sum(index_set),NN);

pre_index=index_set;

while epsilon>GAMMA
    S_set=logical(ones(1,NN));
    candidate_set=S_set;
    candidate_set(index_set)=0;
    
    
   
    while sum(candidate_set)>0
        [kk NN]=size(alpha_zero);
        if length(index_series(index_set))~=kk && sum(pre_index)~=0

            pre_series=index_series(pre_index);
            now_index_set=index_series(index_set);
            [tmp_srt,tmp_idx]=sort([pre_series,setdiff(now_index_set,pre_series)]);
            
            tmp_al=alpha_zero;
            alpha_zero=zeros(length( now_index_set),NN);
            num_new=length(setdiff(now_index_set,pre_series));

            alpha_zero(1:(end-num_new),:)=tmp_al;
            alpha_zero=alpha_zero(tmp_idx,:);
            
            
            pre_index=index_set;
            
        end
        
        
      
        
        
        current_index=index_series(candidate_set);
        current_size=length(current_index);
        rnd_ind=unidrnd(current_size,1,1); 
        this_index=current_index(rnd_ind);
        this_data=matrix_A(index_set,:)';
        p=matrix_A(this_index,:)';
        [inorout,p_prime,alpha_coe,gap]=Spherical_TA(this_data,epsilon,alpha_zero(:,this_index),p);
        if inorout==1
            candidate_set(this_index)=0;

            alpha_zero(:,this_index)=alpha_coe';
        else



            direction=p_prime-p;
            S_index=index_series(candidate_set);
            S_data=matrix_A(S_index,:);

            projected_val=S_data*direction;
            index_2=min(find(projected_val==min(projected_val)));
            index_candidate_2=S_index(index_2);

            
            if index_set(index_candidate_2)==0
                index_set(index_candidate_2)=1;

                candidate_set(index_candidate_2)=0;
            end

        end
    end
%     epsilon
    epsilon=epsilon/2;
end
index=index_series(index_set);

