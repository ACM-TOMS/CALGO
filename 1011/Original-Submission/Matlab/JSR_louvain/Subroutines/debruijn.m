function graph = debruijn(nLabels,dimension)
% DEBRUIJN generates a Debruijn graph with nLabels symbols for a specified 
% dimension (default 1). 
%
%   graph = DEBRUIJN(nLabels)
%     Input nLabels is a positive integer. Returns a DeBruijn graph of 
%     dimension 1 with nLabels symbols.
%
%   graph = DEBRUIJN(nLabels,dimension) 
%     Input dimension is a non-negative integer. Sets the dimension of the 
%     DeBruijn graph. 
%   _______________________________________________________________________
%
%   Description of the DeBruijn graph (see also [1]):
%
%   The well-known DeBruijn graphs are a path-complete graphs : all nodes are
%   possible words of length [dimension] in an alphabet with [nLabels] different
%   letters. There is an edge between node i and node j iff the word of the
%   node i is {a v1 v2 v3 ... v(n-1) } and the word of the node j is 
%   {v1 v2 v3 ... v(n-1) b}, with n the length of the word and a, b, vi are
%    letters of the alphabet. In this case, the label of the edge is 'b'.
%   The number of edges is equal to [nLabels]^([dimension]+1) and the number 
%   of nodes is equal to [nLabels]^[dimension].
%   _______________________________________________________________________
%
%  Thanks to Matthew Philippe for his help for this function.
%
%  [1] J. L. Gross and J. Yellen, Handbook of graph theory (Discrete
%  Mathematics and Its Applications), CRC Press, Boca Raton, FL, 2003.
%
% See also GENERATE_GRAPH

if(nargin < 2 )
    dimension = 1;
end

if(dimension < 0)
    error('Parameter dimension must be positive.')
end

if(not(round(dimension)==dimension))
    error('Parameter dimension must be an integer.')
end

if(not(round(nLabels)==nLabels))
    error('Parameter nLabels must be an integer.')
end

if( nLabels <1 )
    error('Parameter nLabels must be greater or equal than one.')
end

n_edges = nLabels^(dimension+1);

% Memory check
memNeeded = 3*n_edges*8; % supposes double (8 Bytes)
memAvailable = 0.85*available_memory;
if(memNeeded > memAvailable)
    error(['Too much memory needed. Number of edges = [', num2str(n_edges), ']. Memory needed : [',num2str(round(memNeeded/1e6)),'] mB. Memory available : ',num2str(round(memAvailable/1e6)) , ' mB.'] )
end

if(dimension > 5)
    warning('The dimension of the debruijn graph is larger than 5. The computation can be slow.')
end

                                
if(dimension == 0) % One node with only self-loop
    graph.nNodes = 1;
    graph.edges = [ones(nLabels,1),ones(nLabels,1),(1:nLabels)'];
else
    
    graph.nNodes = nLabels^dimension;
    graph.edges = zeros(n_edges, 3); % The edge is [i,j,k] if edge from i
                                     % to j with label k.
    edgeIdx = 0; % Index of the edge.
    for originNode = 1 : graph.nNodes
        originWord = dec2base(originNode - 1, nLabels, dimension); % The node where we come from
        for label = 1 : nLabels
            edgeIdx = edgeIdx +1;
            % All possible words with the same (n-1) first letters.
            destWord = [originWord(2:end), dec2base(label-1,nLabels)]; 
            % We convert the word above into an index -> this is the destination node.
            destNode = base2dec(destWord, nLabels )+1; 
            % We add the edge.
            graph.edges(edgeIdx,:) = [originNode, destNode, label];
        end
    end
end