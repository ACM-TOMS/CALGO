function VV = grVerCover(G);
% l = grVerCover(G);
% Computes all local minimal vertex-covers of a graph
% This algorithm should be fast, if the graph is highly connected.
% If the graph is badly connected, this algorithm is very very slow.
%
% Input:
%   G       matlab-digraph object
%
% Output:
%   l       cell array of indices of vertex covers
%
% Info:
%   Input must be a di-graph object due to technical reasons.
%
% E.g.: G=[1 2; 2 3; 3 4; 2 5];
%       G=digraph(G(:,1),G(:,2));
%       plot(G)
%       vdisp(grVerCover(G))
%
% See also: grCenter
%
% Written by tommsch, 2018

%#ok<*ALIGN>

    VV = grVerCover_new(G); %processes graph componentwise. Should be much faster
end

function VV = grVerCover_new(G);
    %process each connected component on its own
    bins=conncomp(G,'Type','weak','OutputForm','cell');
    VV=cell(1,size(bins,2));
    for i=1:size(bins,2);
        VV{i}=grVerCover_preworker_new(subgraph(G,bins{i}));
        for j=1:size(VV{i},2)
            VV{i}{j}=subsm(VV{i}{j},bins{i},1:numel(bins{i}));
        end
    end


    VV(cellfun(@isempty, VV))=[];%delete empty cells

    VV=mixvector(VV); %make combinations

    VV=arrayfun(@(x) sort(horzcat(VV{:,x})), 1:size(VV,2),'UniformOutput',0); %concatenate combinations

end

function VV = grVerCover_preworker_new(G);

    %plot(G)
    idx=find(any(diff(G.Edges.EndNodes,1,2)==0,2)).'; %indices of self loops
    val=G.Edges.EndNodes(idx,1).';
    G=rmedge(G,idx);
    I=incidence(G); %since self loops are forbidden, make a workaround for it
    for i=val;
        I(i,end+1)=2; %#ok<AGROW>
    end
    %I=full(I); %DEBUG

    %tic
    [VV] = grVerCover_worker_new(I,zeros(1,0),zeros(1,0),zeros(1,0),{}); %make a depth first search, to find all locally minimal vertex covers
    %toc

    VV=uniquecell(VV);

    %these covers can be further simplified, 
    sze=cellfun(@length,VV);
    minlength=min(sze);
    len=cellfun(@(x) length(x(:)), VV);
    for idx=1:size(VV,2);
        if(sze(idx)==minlength); 
            continue; end;
        for j=1:size(VV{idx},2);
            Vtest=VV{idx}([1:j-1 j+1:end]);
            if(searchincellarray(Vtest,VV,1,len));
                VV{idx}=[];
                break;
            end
        end
    end

    VV = VV(~cellfun(@isempty, VV));

end

function [ VV] = grVerCover_worker_new( I, Ec, Eh, V, VV);
    % Ec        Edge_chosen
    % Eh        Edge_hit
    % V         Vertices
    % I         incidence matrix of G

    %choose next edge to process

    Enum=size(I,2); %number of edges
    e=find(~any((1:Enum)'==[Ec Eh],2)); 

    if(isempty(e));  
        if(isequal(1:Enum,Eh));
            VV{end+1}=unique(V); end;
        return; end;
    e=e(1);

    %adjancent vertices of chosen edge
    v=find(I(:,e));

    %call grVerCover_worker three times. 

    Ece=[Ec e]; %union(Ec,e) is slower
    [VV]=grVerCover_worker_new(I, Ece, union(Eh,find(I(v(1),:))), [V,v(1)], VV); %union(V,v(1)) is slower
    if(numel(v)>=2);     
    [VV]=grVerCover_worker_new(I, Ece, union(Eh,find(I(v(2),:))), [V,v(2)], VV); 
    end;
    [VV]=grVerCover_worker_new(I, Ece, Eh                       , V,             VV);
end

function B = subsm(A, newval, oldval)
    B=A;
    for k=1:numel(newval);
        B(A == oldval(k)) = newval(k);
    end;
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   


