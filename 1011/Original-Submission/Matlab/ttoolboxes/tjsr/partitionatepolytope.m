function [chain, MST]=partitionatepolytope(VV, n, verbose);
% [ chain, MST ] = partitionatepolytope( VV, n, verbose );
% Partitions points in $\RR^s$ into clusters of nearby points.
% Also returns a list of consecutive edges. This list is partially ordered.
%
% XX This function needs some more attention. If this function splits up the polytope better, many cores would be used more efficient.
%
% Input:
%   VV          the vertices 
%   n           vertices are split into 2^n parts
%                   The size of each part should be the same, but this does not work
%                   Thus it is good to split it into 2 or 4 times the number of available threads.
%   verbose     verbose level. Default=0.
%
% Output:
%   chain       lists of consecutive edges
%   MST         the graph corresponding to the chain and VV
%
% E.g.: [chain, MST]=partitionatepolytope(normalizematrix(randn(3,1000),'colnorm',4),2,1); 
% 
% Written by: tommsch, 2018

%Test variables
%dim=2; N=200;
%pts=randn(dim,N);
%pts=normalizematrix(pts,'colnorm',2);
%pts=pts+randn(dim,N)*.1;
%W=pdist2(pts.',pts.').';

if(nargin==2); verbose=0; end;


%W=normalizematrix(VV,'colnorm',2);  %Not good: normalize points prior to computing the minspantree
W=VV;

% compute minimal spanning tree
try
    %error; %Debug command to enter catch
    W=exp(pdist2(W.',W.').');
    MST{1}=minspantree(graph(W,'upper'),'Method','dense'); 
catch
    %XX Make numcores many graphs. Right now I only create one graph, thus I have no multithreading anymore.
    nVV=size(VV,2);
    chain{1}=[0:nVV-1; 1:nVV]';
    MST{1}=graph(chain{1}(2:end,1),chain{1}(2:end,2));
    MST{1}.Edges.Weigth=ones(nVV-1,1);
    MST{1}.Nodes.co=VV.';
    MST{1}.Nodes.idx=[1:nVV].';
    vprintf('Out of memory while computing minimum spanning tree. MSP will not be computed.\n','cpr','err');    
    return;
end

clear W;
MST{1}.Nodes.idx=uint32((1:MST{1}.numnodes).'); %save original indices, since matlab reassigns node-ids when making subgraphs
MST{1}.Nodes.co=VV.';
%     plot(MST{1}); %DEBUG


% partitionate graph in subproblems
for i=1:n
    cuts=zeros(0,2);
    for j=1:length(MST);
        center=grCenter(MST{j}.Edges.EndNodes,'edge'); 
        if(isempty(center)); 
            cuts(j,1:2)=[0 1]; 
        else; 
            cuts(j,1:2)=center; 
        end;
    end
    
    for j=1:length(MST);
        if(cuts(j,1)~=0); MST{j}=rmedge(MST{j},cuts(j,1),cuts(j,2)); end;    %#ok<AGROW>
        bins=conncomp(MST{j});
        for k=1:max(bins)
            MST{end+1}=subgraph(MST{j},find(bins==k)); %#ok<AGROW>
        end
        if(max(bins)>2); error('Problem'); end;
        MST{j}=[]; %#ok<AGROW>
    end
    MST = MST(~cellfun(@isempty, MST));
end
clear cuts bins center;

if(verbose>=1);
    %plot graph and paritions
    clf; hold on;
    h=cell(size(MST));
    OCC=zeros(size(MST));
    for j=1:length(MST)

         h{j}=plot(MST{j},'k','LineWidth',2,'Layout','layered','NodeColor',rand(1,3),'MarkerSize',5);

        %h{j}=plot(MST{j},'k','LineWidth',2,'NodeLabel',[],'Layout','layered','NodeColor',rand(1,3),'MarkerSize',5);
        h{j}.XData=VV(1,MST{j}.Nodes.idx);
        h{j}.YData=VV(2,MST{j}.Nodes.idx);
        if(size(VV,1)>=3); h{j}.ZData=VV(3,MST{j}.Nodes.idx); end;
        OCC(j)=MST{j}.numnodes;
    end
    AXIS=axis;
    OCC=sort(OCC);
    plot(linspace(AXIS(1),AXIS(2),length(OCC)),OCC/max(OCC)*AXIS(4));
    axis equal
    vprintf('Size of subtrees: %v\n',OCC);
end
  
%transform vertices of MST into a consecutive list of pairs
chain=cell(size(MST));
for j=1:length(MST)
    chain{j} = dfsearch(MST{j},1,'edgetonew');
    if(~isempty(chain{j}));
        chain{j}=[0 chain{j}(1,1); chain{j}];
    else
        chain{j}=[0 1];
    end

end

end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 
