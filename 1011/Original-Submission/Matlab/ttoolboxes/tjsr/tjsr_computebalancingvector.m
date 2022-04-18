function [vec,alpha]=tjsr_computebalancingvector(type)
% [vec,alpha] = tjsr_computebalancingvector(type)
% This function belongs to tjsr
% Input
%   type        type struct from tjsr, containing
%                       type.counter.numordering                number of orderings
%                       type.counter.numcandidate               number of candidates
%                       type.counter.numnearlycandidate         number of nearly-candidates
%                       type.balancing{i}.cyclictree.V{1}       vertices of the polytope
%                       type.balancing{i}.cyclictree.v0s{1}     leading dual eigenvector
% Output
%   vec        balancing vector.
%   alpha      indices if balancing succeeded.                
%                   alpha > 1      the balancing succeeded.
%                   alpha <=1      balancing failed
%
% Written by tommsch,2018

numo=type.counter.numordering;
numc=type.counter.numcandidate;
numnc=type.counter.numnearlycandidate;

    qij=zeros(numo,numc);
    for i=1:numo; 
        V_i=type.balancing{i}.cyclictree.V{1};
        V_i=num2cell(V_i',2); %transpose, so we dont need to transpose later
        for j=1:numc; 
            vs_j=type.balancing{j}.cyclictree.v0s{1}; %there is only one ordering, since we are balancing. thus: {1} 
            val=cellfun(@(x) abs(x*vs_j), V_i); %balance only against v0s !!!
            qij(i,j)=max(val);
        end; 
    end;
    bij=-log(qij);
    

    %we use the notation from Protasov 2016, p28

    y=zeros(1,numc+1); y(1)=-1; %objective function

    [A,b]=y0yiyj_matrix(numc,bij);
    lb=-inf(1,numc+1); lb(1+numc+1:end)=-inf;
    ub=inf(1,numc+1); ub(1+numc+1:end)=0;
    
    idx=any(~isfinite([A b]),2);
    if(anym(idx));
        type.info.infotext=vprintf('Nan/Inf occured during computation of balancing vector.','cpr','err','imp',[1,type.opt.verbose],'str',type.info.infotext); 
        alpha=0;
        vec=ones(1,type.counter.numordering+1);
        return;
    end;

    if(numc==1);
        vec=1;
        alpha=1;
    else
        if verLessThan('matlab','8.1'); opts=optimset('Display','off','LargeScale','off', 'Simplex','on');  %#ok<LINPROG>
        elseif verLessThan('matlab','8.4'); opts=optimset('Display','off','Algorithm','simplex');           %#ok<LINPROG>
        else; opts = optimoptions(@linprog,'Display','off','Algorithm','dual-simplex');                     %after 2014b
        end        

        %[balancingvector,~,~]=linprog(y,A,b,[],[],[],[],[],opts);    
        [balancingvector,~,~]=linprog(y,A,b,[],[],lb,ub,[],opts);    
        balancingvector=exp(balancingvector);

        alpha=balancingvector(1);
        vec=balancingvector(2:end);
    end

vec(numc+1:numc+numnc,1)=1./max(qij(numc+1:numc+numnc,:),[],2).*(type.cyclictree.orho(numc+1:numc+numnc)).'.*.99; %add factors for extra-vertices
vec(numc+numnc+1:numo,1)=min(1./max(qij(numc+numnc+1:end,:),[],2),1)*type.opt.autoextravertex; %add factors for extra-vertices
    
end

function [A,B] = y0yiyj_matrix(maxj, bij)
    L=maxj+1;
    A=zeros(maxj^2+4*maxj,L);
    B=zeros(maxj^2+4*maxj,1);
    rowidx=1;
    for i=2:maxj+1
        for j=2:maxj+1
            if(i==j); continue; end;
            A(rowidx,1:L)=zeros(1,L); 
            A(rowidx,1)=1;
            A(rowidx,i)=1;
            A(rowidx,j)=-1;
            B(rowidx,1)=bij(i-1,j-1);
            rowidx=rowidx+1;
        end
    end
    A(rowidx:end,:)=[];
    B(rowidx:end,:)=[];
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   