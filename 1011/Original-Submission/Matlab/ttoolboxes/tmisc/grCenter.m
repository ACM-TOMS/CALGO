function [center, weight]=grCenter(varargin)
% [center, weight]=grCenter( E, [options] )
% Finds the center of a tree
% Iteratively deletes leaves, until the center is found.
%
% Input:
%   E       List of edges (N x 2)-array
%
% Options:
%   'edge'      counts vertices to return an edge in all cases. 
%               This algorithm does not work for trees with multiple components
%
% Output:
%   center      index of vertices/pairs of vertices of the edge which are the center
%               If tree has more than one component, center is a cell array
%   weight      (behaviour may be changed) distance to mostouter leave
%
% Note: 
%   If E is not a tree, the behaviour is undefined
%
% E.g.: G=[1 2;2 5; 5 6; 2 3; 3 4; 4 7; 4 8; 8 9]; 
%       plot(graph(G(:,1),G(:,2))); 
%       grCenter(G)
%       grCenter(G,'edge')
%
% Written by: tommsch, 2018

% XX Abort loop if graph is not a tree

%#ok<*ALIGN>

E=varargin{1}; varargin(1)=[];
if(size(E,1)==1);  
    center=E;  
    weight=[0 0]; 
    return; end;
if(size(E,1)==0); 
    center=[]; 
    weight=[]; 
    return; end;

[edgeflag,varargin]=parsem('edge',varargin);
if(edgeflag); 
    [center, weight] = grCenter_edge(E,varargin{:});
else
    [center, weight] = grCenter_vertex(E,varargin{:});
end

end

function [center, weight]=grCenter_edge(E,varargin)

%G=graph(E(:,1),E(:,2));
V=unique(E(:)); %list of vertices
%num=numel(V);
weight=zeros(size(V));
while(~isempty(E))
    
    %identify leaves
    [occ,idx]=hist(E(:),V); %occurences, vertex-idx
    idx=idx(find(occ==1)); %#ok<FNDSB>
    weight(idx)=weight(idx)+1;
    
    %plot(G,'NodeLabel',weight)
    
    %sum up leaves
    
    edges=arrayfun(@(x) find(any(E==x,2)),idx,'UniformOutput',0);
    edges=cell2mat(edges);
    E2=E(edges,:);
    for i=1:size(edges,1)
        val=weight(E2(i,1))+weight(E2(i,2));
        
        if(E2(i,1)==idx(i));
            weight(E2(i,2))=val;
        else
             weight(E2(i,1))=val;
        end

        %plot(G,'NodeLabel',weight)
    end
    
    
    %remove edges
    %E(edges,:)
    E(edges,:)=[];
end


[~,center]=sort(weight,'descend');
center=center(1:2);

%h=plot(G,'NodeLabel',weight);
%highlight(h,center);
end

function [center, weight]=grCenter_vertex(E,varargin)


weight=zeros(size(E,1),1);
%plot(graph(E(:,1),E(:,2)))

k=1; %counter for weight
kidx=1;

V=unique(E(:)); %list of vertices
center={};
possiblecenter=[];
while(true);
    [occ,idx]=hist(E(:),V); %occurences, vertex-idx
    idx=idx(find(occ==1)); %#ok<FNDSB>
    
    val=possiblecenter(occ(possiblecenter)==0);
    if(~isempty(val)); 
        center=[center num2cell(val.')]; %#ok<AGROW>
    end;
    
    if(isempty(E)); 
        break; end;
    possiblecenter=setdiff(V(occ~=0),idx);
    edges=arrayfun(@(x) find(any(E==x,2)),idx,'UniformOutput',0);
    
    edges=cell2mat(edges);
    
    %test for multiple edges ==> meaning: this edge is in the center
    [~, I] = unique(edges, 'first');
    val = 1:size(edges,1);
    val(I) = [];
    val=edges(val,:);
    if(~isempty(val));
        center{end+1}=E(val,:); %#ok<AGROW>
    end
    
    
    weight(kidx:kidx+size(edges,1)-1)=repmat(k,size(edges,1),1);
    
    kidx=kidx+size(edges,1);
    k=k+1;
    E(edges,:)=[];
    
end

if(size(center,2)==1); center=center{1}; end;


end