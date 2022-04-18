function type = tjsr_generatenewvertex(M, type);
% type = tjsr_generatenewvertex(M, type);
% This function belongs to tjsr!
% Generates all possible new vertices
%
% XX Implement that testeigenplane can find new s.m.p.-candidate
%
% Written by tommsch, 2018


for i=1:type.counter.numordering %iterate over all trees %i is the tree
    
    %get indices of vertices which get children ad construct orderings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parentidx = ~isnan(type.cyclictree.norm{i});                                %get vertices which can become parent
    parentidx = parentidx & ~(type.cyclictree.norm{i}<1-type.opt.epspolytope);  %get vertices which are outside
    parentidx = parentidx & ~type.cyclictree.status{i}==1;                      %test if childrean are already constructed
    %construct indices of parents
    parent=1:length(parentidx); 
    parent=parent(parentidx);
    parent=repmat(parent,1,type.counter.nummatrix);
    %make all new orderings
    oo=type.cyclictree.o{i}(:,parentidx);  %oo: all possible new orderings
    noo=size(oo,2); 
    oo=repmat(oo,1,type.counter.nummatrix);
    oo=[oo; reshape(repmat(1:type.counter.nummatrix,noo,1),1,[])]; %#ok<AGROW>
    oo=removezero_compact_fast(oo);
    
    %find orderings which are already in the tree somewhere
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(type.cyclictree.smpflag(i)~=2 && ~isempty(oo));
        oclass=type.cyclictree.oclass{i};
        if(size(oo,1)>size(oclass,1)); 
            oclass(size(oo,1),end)=0; % make oo and oclass the same length
        elseif(size(oo,1)<size(oclass,1)); 
            oo(size(oclass,1),end)=0; 
        end;
        removeidx=ismember(oo',oclass','rows').'; %indices of orderings are in oclass
        oo(:,removeidx)=[]; %delete those
    else
        removeidx=false(1,size(oo,2));
    end

    %set all properties for the new vertices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    szeoo=size(oo,2);
    type.cyclictree.status{i}(parent)=1;  %set status of all parents to one (i.e. they are parents)
    parent=parent(~removeidx);
    normparent=repmat(type.cyclictree.norm{i}(parentidx),1,type.counter.nummatrix);
    normparent=normparent(~removeidx);
    type.cyclictree.level{i}(end+1:end+szeoo)=length(type.cyclictree.L{i})+1;
    type.cyclictree.status{i}(end+1:end+szeoo) =        zeros(1,szeoo); %set status for all new entries to zero (i.e. they are children)
    type.cyclictree.normparent{i}(end+1:end+szeoo) =   normparent;
    type.cyclictree.parent{i}(end+1:end+szeoo)=parent;

    %the following variables are needed due to parallelization
    [vv,vs]=deal(zeros(type.info.dim,szeoo));   %local variable where coordinates ofnew vertices are saved
    rho_val=zeros(1,szeoo);                     %local variable where spectral radii of matrix products are saved
    vsflag=type.opt.testeigenplane>-inf || type.opt.balancing;  %if this flag is set, Vs is computed
    testspectralradius=type.opt.testspectralradius;             %if this flag is set, spectral radii are tested
    v0=type.cyclictree.v0{i};                   %starting vector of cyclic root
    v0s=type.cyclictree.v0s{i};                 %dual starting vector of cyclic root
    parfor j=1:szeoo
        newM=tbuildproduct_fast(M, oo(:,j));
        vv(:,j)=newM*v0;
        if(vsflag); vs(:,j)=newM*v0s; end;
        if(testspectralradius); rho_val(j)=trho(newM)^(1/nnz(oo(:,j))); end;
    end    
    
    type.cyclictree.normest{i}(end+1:end+szeoo)=NaN; %norm-estimate is done in another function, because I want to be able to estimate all norms in each round again.
    type.cyclictree.norm{i}(end+1:end+szeoo)=NaN;
    if(testspectralradius); type.cyclictree.rho{i}(end+1:end+szeoo)=rho_val; end;
    if(vsflag); type.cyclictree.Vs{i}(:,end+1:end+szeoo)=vs; end;
    type.cyclictree.V{i}(:,end+1:end+szeoo)=vv;
    type.cyclictree.L{i}(end+1)=szeoo;
    type.cyclictree.o{i}(1:size(oo,1),end+1:end+szeoo)=oo;
    
    %rholeqval: remove all vertices whos matrix products have spectral radius greater equal 1+10*eps %DEBUG CODE - Can be removed safely
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rholeqval=type.opt.rholeqval;
    if(any(rho_val) && rholeqval)
        idx=type.cyclictree.rho{i}>=rholeqval;
        if(any(idx))
            type.cyclictree.V{i}(:,idx)=0;
            type.cyclictree.norm{i}(:,idx)=-.5; %set to a value which can be recognized
            type.cyclictree.normest{i}(:,idx)=0;
            rho_val=0;
            type.info.infotext=vprintf('Removed %i matrices due to ''rholeqval''.\n',nnz(idx),'cpr','err','imp',[1,type.opt.verbose],'str',type.info.infotext);
        end
    end
    
    %test stopping-criterion: eigenplane 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(type.opt.testeigenplane>-inf && type.cyclictree.smpflag(i)==0);
        err_eigenplane=zeros(1,szeoo);
        parfor j=1:size(oo,2)  
            err_eigenplane(j)=abs( dot(v0s,vv(:,j)) ); %'abs'/'real' are necessary for the case real/complex - eigenvectors %COMPLEX CASE
        end
        idx_eigenplane=find(err_eigenplane>1-type.opt.testeigenplane,1);
        if(any(idx_eigenplane))
            type.info.errorcode=110; % 110 = TJSR_TESTEIGENPLANE;
            type.info.errorinformation={max(err_eigenplane)};
        end
    end    
    
    %test stopping-criterion: spectral radius
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(type.opt.testspectralradius>-inf);
        [val,idx]=max(rho_val);
        if(val>type.opt.delta*(1-type.opt.testspectralradius));
            oo_rho=reducelength(oo(:,idx));
            loo_rho=length(oo_rho);
            type.info.errorcode=60; % 60 = TJSR_BETTERORDERINGFOUND;
            type.info.errorinformation={trho(tbuildproduct(type.M_original,oo_rho))^(1/loo_rho),oo_rho }; %{rho, ordering}
            %for j=1:size(type.cyclictree.ordering,2)
            %    if(loo_rho>max(cellfun(@(x) size(x,1), type.cyclictree.ordering)) || ...
            %        ismember(...
            %        [oo_rho; zeros(size(type.cyclictree.ordering{j},1)-loo_rho,1)]',...
            %        type.cyclictree.ordering{j}','rows')...
            %      )
            %    end
            %end
        end
    end
end



end

function c=removezero_compact_fast(c)
    %remove all zero rows
    idx=~any(c,2); 
    c(idx,:)=[];
    
    %compact columns
    len=size(c,1);    
    for i=1:size(c,2);
        idx=c(:,i)~=0;
        c(:,i)=[c(idx,i); zeros(len-nnz(idx),1)];
    end
    
    %remove all zero rows again
    idx=~any(c,2); 
    c(idx,:)=[];
end




function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   