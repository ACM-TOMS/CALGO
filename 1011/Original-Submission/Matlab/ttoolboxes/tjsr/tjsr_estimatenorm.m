function type = tjsr_estimatenorm( type );
% This function belongs to tjsr
% Written by: tommsch, 2018
%
% XX This function should be included in computepolytopenorm



    if(type.opt.fastnorm<1); 
        return; end; %do not estimate if option is set

    %check if we test old vertices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(type.opt.testoldvertex==0); %we do not test old vertices
        selectoldvertexflag=0;
    elseif(type.opt.testoldvertex==1) %if number of vertices increases much, we test old vertices
        REESTIMATE=2; %by what factor the number of vertices of the polytope must grow, until we reestimate all vertices
        MINNUM=500; %minimum number of old vertices to test. If set too low, in the first runs old vertices are tested always
        val=log2(sum(cumsum(vertcat(type.cyclictree.L{:}),2),1))/log2(REESTIMATE);
        [~,~,val]=unique(ceil(val));
        val=diff(val');
        val=val(end);
        if(val && sum([type.cyclictree.L{:}])>MINNUM); 
            selectoldvertexflag=1; 
        else; 
            selectoldvertexflag=0; end;
    elseif(type.opt.testoldvertex==2) %we always test old vertices
        selectoldvertexflag=1;
    else
        error('Wrong value for ''testoldvertex''.'); end;
    if(selectoldvertexflag)
        type.info.infotext=vprintf('Estimate old vertices. ','imp',[2,type.opt.verbose],'str',type.info.infotext); end;


    %Estimate norms
    %%%%%%%%%%%%%%%%%%%%%%%%%
    VV = tjsr_getpolytope(type); %get polytope for which we compute the norms    
    incounter=0; %number of points which are prooven to be inside
    outcounter=0; %number of points which are prooven to be outside
    maxlevel=size(type.cyclictree.L{1},2);
    nestimate=0; %number of estimated vertices in total
    
    for i=1:type.counter.numordering %iterate through all trees
        %select vertices which we want to estimate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx_all = type.cyclictree.level{i}==maxlevel; %select all new vertices
        if(selectoldvertexflag); %also select old vertices if there are enough new vertices in the polytope
            OFFSET = 2; %Heuristic value
            %only test vertices which are added long ago (OFFSET many iterations before)
            %only test vertices which are not inside (this includes untested vertices)
            idx = type.cyclictree.level{i} < type.cyclictree.level{i}(end)-OFFSET  &  type.cyclictree.norm{i} >= 1-type.opt.epsequal;
            idx_all = idx_all | idx; end;

        %Estimate norm
        %%%%%%%%%%%%%%%%%
        pts=type.cyclictree.V{i}(:,idx_all); %vertices to test
        nestimate=nestimate+size(pts,2); %increase counter
        type.info.infotext=vprintf('Estimate %i vertices for tree %i.\n',size(pts,2),i, 'imp',[4,type.opt.verbose],'str',type.info.infotext); 
        [normest,region]=estimatepolytopenorm(pts,VV,[],type.info.algorithm,type.opt.epsequal);

        %Save estimate
        %%%%%%%%%%%%%%%%%
        type.cyclictree.normest{i}(idx_all)=normest;
        idx_all=find(idx_all);

        if(type.opt.fastnorm>=1) %save norms for points prooven to be inside
            idx=normest<region(1);
            idx=idx_all(idx);
            type.cyclictree.norm{i}(idx)=TJSR_INSIDE;
            incounter=incounter+length(idx); end;

        if(type.opt.fastnorm>=2);  %save norms for points prooven to be outside
            idx=normest>region(2);
            idx=idx_all(idx);
            type.cyclictree.norm{i}(idx)=min(type.cyclictree.norm{i}(idx),TJSR_OUTSIDE); 
            outcounter=outcounter+length(idx); end;
    end;
    
    %Text Output
    %%%%%%%%%%%%%%%%%%%%%%%%
    if(nestimate);  
        type.info.infotext=vprintf('Estimate %i vertices. ', nestimate, 'imp', [3,type.opt.verbose], 'str', type.info.infotext);    end
    if(incounter);  
        type.info.infotext=vprintf('In: %i, ', incounter ,'imp', [2,type.opt.verbose], 'str', type.info.infotext);                  end;
    if(outcounter); 
        type.info.infotext=vprintf('Out: %i, ', outcounter,'imp', [2,type.opt.verbose], 'str', type.info.infotext);                 end;
    if(nestimate);  
        type.info.infotext=vprintf('\n', 'imp', [3,type.opt.verbose], 'str', type.info.infotext);                                   end
end

function val = TJSR_INSIDE;                         val=.5; end %Arbitrary value smaller than 1-epsilon
function val = TJSR_OUTSIDE;                        val=inf;  end %Arbitrary value larger than 1


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   