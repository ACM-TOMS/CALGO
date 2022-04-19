function [type]=tjsr_makecyclicroot(M, type);
% type = tjsr_makecyclicroot(M, type)
% This function belongs to tjsr!
%
% Constructs .cyclictree as it shall be, after the first iteration
% Make sets V, o, L, Vs, Rho
% V contains the vertices of the polytope
% o(:,i) contains the candidates of the matrices how to obain the vertex V(:,i)
% k counts the number of processed steps
% L(k) contains the number of vertices added in the k^th step
% V and o is such that: buildProduct(A,o(:,n))*V(:,1)=V(:,n)
% vs are the leading eigenvectors of PI_j
%
% Written by tommsch, 2018
    
    %%%%Good matrices to test this function:
    % tjsr(tgallery('gallery_randcolu_1',10,20,7),'invariantsubspace','none','verbose',3,'nearlycandidate',0.99,'maxsmpdepth',10,'nobalancing','extravertex',{[1 1 1 1 1 1 1 1 1 1]'})
    % tjsr(tgallery('rand_rho',2,2,7),'invariantsubspace','none','verbose',3,'nearlycandidate',0.99,'maxsmpdepth',10,'nobalancing','extravertex',{[1;1]})
    %
    
    %Entries for each Level (Iteration)
    type.cyclictree.L=cell(int32(1),type.counter.numordering);  [type.cyclictree.L{:}]=deal(int32(zeros(1,0)));
    type.cyclictree.livingvertex=zeros(1,0);
    type.cyclictree.timelvl=zeros(1,0); 
    type.cyclictree.normlvl=zeros(1,0); 
    %type.cyclictree.flatness=zeros(1,0);
    
    %Entries for vertices
    type.cyclictree.level=      cell(1,type.counter.numordering); %depth of the vertex == distance of vertex to cyclic root    
    type.cyclictree.norm=       cell(1,type.counter.numordering); %norm of the vertex corresponding to the polytope (generated so far).     
    type.cyclictree.normest=   cell(1,type.counter.numordering); %estimate of the polytope norm based on crude computations
    type.cyclictree.normparent=cell(1,type.counter.numordering); %norm of the parent vertex
    if(type.opt.maxtreedepth>500);  type.cyclictree.o=repcell(int32(zeros(20,1)),1,type.counter.numordering); %the orderings of the vertices
    else;                           type.cyclictree.o=      repcell(int32(zeros(type.opt.maxtreedepth,1)),1,type.counter.numordering); 
    end
    type.cyclictree.parent=     cell(1,type.counter.numordering); %index of parent vertex
    type.cyclictree.rho=        cell(1,type.counter.numordering); %spectral radius of the matrix product corresponding to the vertex
    type.cyclictree.status=     cell(1,type.counter.numordering); %flag if vertex has children or not
    type.cyclictree.V=          repcell(zeros(type.info.dim,1),1,type.counter.numordering);  %coordinates of the vertices
    type.cyclictree.Vs=         repcell(zeros(type.info.dim,1),1,type.counter.numordering);  %coordinates of the *vertices
    
    for i=1:type.counter.numordering
        if(type.cyclictree.smpflag(i)~=2) %for candidates and nearlycandidates
            norderingi=size(type.cyclictree.ordering{i},2);
            for k=1:norderingi
                orderingk=removezero(type.cyclictree.ordering{i}(:,k),'all'); %the ordering which gets processed now
                szek=size(orderingk,1);
                type.cyclictree.o{i}(1,1)=int32(0); %contains the ordering for the vertex 
                type.cyclictree.V{i}(:,1)=type.cyclictree.v0{i}; %save zeroth entry
                type.cyclictree.Vs{i}(:,1)=type.cyclictree.v0s{i};
                for j=1:szek-1 %I start from j=1, since j=0 is the empty product which is already added before. I only must add it once.
                    newo=orderingk(1:j);
                    newv=makepositive(tbuildproduct(M,newo)*type.cyclictree.v0{i});
                    newvs=makepositive(tbuildproduct(M,newo)'*type.cyclictree.v0s{i});
                    type.cyclictree.o{i}(1:j,end+1)=int32(newo); %contains the ordering for the vertex 
                    type.cyclictree.V{i}(:,end+1)=newv; 
                    type.cyclictree.Vs{i}(:,end+1)=newvs;
                    if(~isempty(newo)); 
                        type.cyclictree.oclass{i}(1:j,end+1)=int32(newo);   %add vertex to oclass 
                    end;
                end
            end
            %remove duplicates
            if(size(type.cyclictree.o{i},2)>=2); 
                [~,idx,~]=unique(type.cyclictree.o{i}','rows');
                type.cyclictree.o{i}=type.cyclictree.o{i}(:,idx);
                type.cyclictree.V{i}=type.cyclictree.V{i}(:,idx);
                type.cyclictree.Vs{i}=type.cyclictree.Vs{i}(:,idx);
            end
            
            type.cyclictree.level{i}(1:size(type.cyclictree.o{i},2))=int32(1);
            type.cyclictree.norm{i}(1:size(type.cyclictree.o{i},2))=inf; %it is inf by convention, since it is the norm of the point wrt. no polytope
            type.cyclictree.normest{i}(1:size(type.cyclictree.o{i},2))=inf; %it is inf by convention, since it is the norm of the point wrt. no polytope
            type.cyclictree.normparent{i}(1:size(type.cyclictree.o{i},2))=inf; %by definition
            type.cyclictree.parent{i}(1:size(type.cyclictree.o{i},2))=int32(0); %vertices in the root have no parents by definition, since otherwise they will not get constructed in the algorithm later
            type.cyclictree.status{i}(1:size(type.cyclictree.o{i},2))=int8(0); %no children by definition, since otherwise they will not get constructed in the algorithm later
            type.cyclictree.rho{i}(1:size(type.cyclictree.o{i},2))=trho(tbuildproduct(M,type.cyclictree.ordering{i}(:,1)))^(1/length(type.cyclictree.ordering{i}(:,1))); %1 for candidates, <1 for nearlycandidates
            type.cyclictree.L{i}(1)=int32(size(type.cyclictree.o{i},2));
        else %for extravertices

            type.cyclictree.level{i}=int32(1);
            type.cyclictree.norm{i}(1)=Inf; %it is inf by convention, since it is the norm of the point wrt. no polytope
            type.cyclictree.normest{i}(1)=inf; %it is inf by convention, since it is the norm of the point wrt. no polytope
            type.cyclictree.normparent{i}(1)=inf; %by definition
            type.cyclictree.parent{i}(1)=int32(0); %no parent
            type.cyclictree.rho{i}(1)=NaN; %it is NaN, since it stems from no matrix product
            type.cyclictree.status{i}(1)=int8(0); %no children
            type.cyclictree.V{i}(:,1)=type.cyclictree.v0{i};
            type.cyclictree.L{i}(1)=int32(1);
        end
    end
    
    %set values for level-things
    type.cyclictree.timelvl(1)=0;
    type.cyclictree.normlvl(1)=inf;
    %type.cyclictree.flatness(1)=flatness([type.cyclictree.V{:}],type.info.algorithm);
    type.cyclictree.livingvertex(1)=size([type.cyclictree.V{:}],2);
    
    %simplifyo oclass
    for i=1:type.counter.numordering
        type.cyclictree.oclass{i}=unique(type.cyclictree.oclass{i}.','rows').'; end;
    
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 