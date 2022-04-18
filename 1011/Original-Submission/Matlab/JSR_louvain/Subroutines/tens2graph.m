function graph = tens2graph(tens)
%
% graph = TENS2GRAPH(tens)
%   Input tens is a 3-D tensor with only binaries values. 
%   Convert a tensor to a graph structure used in the JSR toolbox. If 
%   tens(i,j,k) == 1, then there exists a edge from node i to node j with 
%   label k in the graph.
%
%   The output has the following structure.
%               1) graph.nNodes     Contains the number of nodes in the graph.
%                                   Numerotation begins at 1.
%               2) graph.edges      nedgesx3 Matrix of integer such that
%                                   graph.nedges(m,:) = [i,j,k] is the mth
%                                   edges of the graph, and the edge comes
%                                   from node i, goes to node j, and have
%                                   the label M{k}.
%                                   Numerotations begin at 1.
%
%   EXAMPLE : Random tensor.
%   >> tens = randi(2,3,3,4)-1; % random binary tensor of size (3,3,4)
%   >> graph = tens2cell(tens);


sizeMatrix = size(tens,1);
if size(tens,2) ~= sizeMatrix
    error('Input matrices must be square');
end
nLabel = size(tens,3);
nEdges = sum(sum(sum(tens)));


memNeeded = nEdges*3*8;
memAvailable = 0.85*available_memory;
if(memNeeded > memAvailable)
    error(['Too much memory needed. Number of edges = [', num2str(n_edges), ']. Memory needed : [',num2str(round(memNeeded/1e6)),'] mB. Memory available : ',num2str(round(memAvailable/1e6)) , ' mB.'] )
end

graph.edges = zeros(nEdges,3);
graph.nNodes = sizeMatrix;

edgeCounter = 0;
for k=1:nLabel
    for i = 1:graph.nNodes
        for j = 1:graph.nNodes
            
            if not(tens(i,j,k) == 0 || tens(i,j,k) == 1)
                error('The tensor must contain only binaries values.')
            end
            
            if tens(i,j,k) == 1
                edgeCounter = edgeCounter+1;
                graph.edges(edgeCounter,:) = [i,j,k];
            end
        end
    end
end