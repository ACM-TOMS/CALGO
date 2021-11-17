function Icomp = complementarityConstraint(maxdeg, factor, clique, randSeed)

% Kojima
    n = size(clique.Set,2);
    rng('default');
    if nargin < 4; randSeed = 2017; end;
    rng(randSeed);
    if maxdeg < 2; maxdeg = 2; end

    blksizes = sum(clique.Set, 2);
    Icomp = false(factor * sum(blksizes), n);
    ptr = 1;
    for c = 1:size(clique.Set, 1)
        blkidx = find(clique.Set(c, :));
        blksize = blksizes(c);
        numCompInBlk = ceil(factor * blksize);
        for k = 1:numCompInBlk
            degComp = 1 + ceil(rand() * (maxdeg - 1));
            if (blksize > 3) && (blksize < 2 * degComp)
                continue
            end
            idx = [];
            while (length(unique(idx)) <= 1) 
                idx = ceil(rand(1, degComp) * blksize);
            end
            Icomp(ptr, blkidx(idx)) = true;
            ptr = ptr + 1;
        end
    end
    Icomp(ptr:end, :) = [];
    Icomp = unique(Icomp, 'rows');
end