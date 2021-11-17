function [triangul, blocks, perm] =  tpermTriangul(varargin)
% [ triangul, blocks, perm ] =  tpermTriangul( M, [options] )
%  Looks for a permutation of the rows and columns of the matrices in M that jointly block-triangularize them.
%
% [TRIANGUL, BLOCKS, PERM] =  tpermTriangul(M)
%    For a cell-array M of nonnegative square matrices. If it exists, 
%    finds a permutation PERM such that :
%                  
%                   blocks{1}{i}    *          *       *    . . .    * 
%                      0    
%                               blocks{2}{i}   *       *             .
%                      0            0                                .
%                      .            .       blocks{3}{i}             .
% M{i}(PERM,PERM) =    .            .          0         .           
%                      .            .                        .  
%                                                               .    *
%                      0            0          0     . . .       blocks{q}{i}  
%
%
% Options:
%   'verbose',val       Verbose level
%                                             
% REFERENCES
%    R.Jungers, 
%      "The Joint Spectral Radius: Theory and Applications" 
%      Vol. 385 section 1.2.2.5 in Lecture Notes in Control and Information
%      Sciences, Springer-Verlag. Berlin Heidelberg, June 2009
%    R.Jungers, V.Protasov, V.Blondel,
%      "Efficient algorithms for deciding the type of growth of products of integer matrices"
%      Linear Algebra and its Applications, 428(10):2296-2311, 2008
% 
% Modified by: tommsch, 2019

M=varargin{1}; varargin(1)=[];
verbose=parsem({'verbose','v'},varargin,1);

n = size(M{1},1);
m = length(M);

Sum = sparse(n,n);

for i=1:m
     Sum = Sum + abs(M{i});
end

Sum = Sum - diag(diag(Sum));

% Adjancy matrix of the graph
G = sparse(isAlways(Sum>0,'Unknown','true'));


% Find the strongly connected components
[S,C] = tgraphSCC(G);

% Visualization feature
%
if(verbose>=2);
    % Mark the nodes for each component with different color
    h = view(biograph(G));
    colors = jet(S);
    for i = 1:numel(h.nodes)
        h.Nodes(i).Color = colors(C(i),:);
    end
end

if S==1
    triangul = 0;
    blocks = M;
    perm = 1:n;
else
    triangul = 1;
    blocks = cell(1,S);
    perm = zeros(1,n);
    ind = 0;

    for i=S:-1:1
        ncon = sum(C==i);
        perm(ind+(1:ncon))=find((C==i));
        for imat = 1:m
            blocks{S-i+1}{imat} = M{imat}(perm(ind+(1:ncon)),perm(ind+(1:ncon)));
        end
        
        ind = ind +ncon;
    end
end

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 