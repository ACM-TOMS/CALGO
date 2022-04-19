function chain = findLongChain(supports, eq0Var, chain)
    [m, n] = size(supports);
    maxdeg = max(sum(supports, 2));
    unused = ~(eq0Var | any(chain, 2));%true(m, 1);
    
    col = zeros(m, 1); % idx of chain
    row = zeros(m, 1); % elements in chain
    ch = 1;
    it = 1;

    function addElement(r)
        col(it) = ch;
        row(it) = r;
        unused(r) = false;
        it = it + 1;
    end

    while any(unused)
        %% starting node of chain
        s = find(unused, 1);
        if sum(supports(s, :)) > ceil(maxdeg/2); break; end
%         if sum(supports(s, :)) >= maxdeg; break; end
        addElement(s);
        
        %% next node
        if true
            while true
                NZVar = logical(supports(s, :));
                ns = sum(NZVar); 
                %allocate boolean matrix in dense format
                tmp = false(m, ns); 
                %substitute into the dense format directly (avoid creating dense matrix in sparse format)
%                 tmp(:) = (supports(:, NZVar) >= supports(s, NZVar)); 
                tmp(:) = (supports(:, NZVar) >= repmat(supports(s, NZVar), size(supports, 1), 1)); 
                t = find(unused & all(tmp, 2), 1);
                if isempty(t); break; end
                addElement(t);
                s = t;
            end
        else
            % supports >= supports(s, :) makes dense matrix in sparse format
            % that is too inefficient.
            t = find(unused & all(supports >= supports(s, :), 2), 1);
            while ~isempty(t)
                addElement(t);
                s = t;
                t = find(unused & all(supports >= supports(s, :), 2), 1);
            end
        end
        ch = ch + 1;
    end
    row(it:end) = [];
    col(it:end) = [];
    chain2 = sparse(row, col, true, m, ch-1);
    chain2 = chain2(:, sum(chain2, 1) ~= 1);
    chain = [chain, chain2];
end