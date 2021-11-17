function [scnt,C] = tgraphSCC(G)
% 
% tgraphSCC  Finds the strongly connected components of graph 
%           G using Tarjan's algorithm.
%
% [SCNT,C] = tgraphSCC(G)
%
%  G must be the adjency matrix of a directed graph, i.e. if entry (i,j) is
%    non-zero there is an edge going from i to j
%
%  SCNT is the number of strongly connected components 
%
%  C is a vector indicating to which component each node belongs
%
%
% REFERENCES: 
%    R.Tarjan, 
%      "Depth first search and linear graph algorithms", 
%      SIAM J. on Computing, 1(2):146-160, 1972 
%    R.Sedgewick, 
%      "Algorithms in Java, Third Edition, Part 5: Graph Algorithms",
%      Section 19.8, Addison-Wesley, 2003
%
% Written by: Jungers, The JSR Toolbox

global S


n = length(G);
Adj = cell(n,1);

for ivert = 1:n
    Adj{ivert} = find(G(ivert,:)>0);
end

scnt=0;
cnt = 0;

S.nel=0;
S.stack = {};
S.isIn = zeros(1,n);

low = -ones(1,n);
number = -ones(1,n);
C = -ones(1,n);

for iv = 1:n
    if (number(iv)==-1)
        strongConnect(iv);
    end
end

    function strongConnect(v)
        cnt = cnt+1;
        low(v)=cnt;
        number(v)=cnt;
        add(v);
        
        for w = Adj{v}
            if (number(w)==-1)
                strongConnect(w);
                low(v) = min(low(v),low(w));
            elseif (number(w)<number(v))
                if (S.isIn(w))
                    low(v) = min(low(v),number(w));
                end
            end
        end
        if (low(v) == number(v))
            cont = 1;
            scnt = scnt+1;
            
            while (S.nel ~= 0 && cont)
                top = S.stack{1};
                if (number(top) >= number(v) )
                    w = pop();
                    C(w) = scnt;
                else 
                    cont = 0;
                end
            end
        end
    end
end

function add(v)
global S
  S.nel = S.nel+1;
  S.stack = [{v} S.stack]; 
  S.isIn(v)=1;
end

function  v = pop()
global S
  if (S.nel~=0)
      v = S.stack{1};
      S.nel = S.nel-1;
      S.stack(1) = [];
      S.isIn(v)=0;
  else
      v = [];
  end
    
end
