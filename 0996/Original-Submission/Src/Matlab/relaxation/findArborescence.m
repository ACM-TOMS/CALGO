function chain = findArborescence(supports, eq0Var, chain)
%%FINDARBORESCENCE Find a spanning arborescence
% NOTE: Currently, arborescence cannot be dealt within BP.m function.
% This function adds new edges into chains.
% The new edges will be dealt as general inequality constraints in 
% interior point methods.
%%
rng(1)
[m, n] = size(supports);
deg = sum(supports, 2);
num = accumarray(full(deg)+1, 1);
S = mat2cell(supports, num, n);
isConnected = eq0Var | any(chain, 2);
for jj = 1:size(chain, 2)
    isConnected(find(chain(:, jj), 1)) = false;
end
isConnected = mat2cell(isConnected, num, 1);
tree = sparse(false(m, m-1));
tree = mat2cell(tree, num, m-1);
%%
treeidx = 1;
for d = (max(deg):-1:1) + 1
    for ii = 1:size(S{d}, 1)
        if isConnected{d}(ii)
            continue
        end
        u = S{d}(ii, :);
        for jj = randperm(size(S{d-1}, 1))
            v = S{d-1}(jj, :);
            if all(u >= v)
                tree{d}(ii, treeidx) = true;
                tree{d-1}(jj, treeidx) = true;
                treeidx = treeidx + 1;
                isConnected{d}(ii) = true;
                break;
            end
        end
    end
end
tree = cell2mat(tree);
tree(:, ~any(tree, 1)) = [];
chain = [chain, tree];
end

