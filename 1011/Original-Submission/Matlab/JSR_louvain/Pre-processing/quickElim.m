function [newblocks, elim, bounds] = quickElim(blocks)
%
% [NEWBLOCKS, ELIM, BOUNDS] = QUICKELIM(BLOCKS)
%      If BLOCKS is a cell array of sets of blocks (for instance from
%      permTriangul or jointTriangul), this method founds bounds by a few
%      bruteForce iterations in order to take out the blocks that will not
%      play a role in jsr(M) = max(jsr(blocks{1},...,blocks{q})).
%
%  NEWBLOCKS is the cell array of the sets of blocks that could not be
%            dismissed
%  ELIM is a vector of length = length(BLOCKS), ELIM(ind)=1 if one could
%       dismiss blocks{ind},
%
%  BOUNDS contains on its rows the bounds found for the corresponding set
%         of BLOCKS
%
% Hence, ELIM(ind)=1 iff there is some set of blocks BLOCKS{j} s.t. 
% bounds(ind,2)<bounds(j,1)
%  

nblo = length(blocks);
lbBlo = zeros(1,nblo);
ubBlo = zeros(1,nblo);

for iblo = 1:nblo
    Msub = blocks{iblo};
    m = length(Msub);
    if (size(Msub{1},1)==1)
        lbBlo(iblo) = max(abs(cell2mat(Msub)));
        ubBlo(iblo) = lbBlo(iblo);
    else
        mdbrute  = ceil(log(10)/log(m));
        optsQuickBrute = jsrsettings('verbose',0,'logfile',0,'bruteforce.maxdepth',mdbrute);
        bQuickBrute = jsr_prod_bruteForce(Msub,optsQuickBrute);
        lbBlo(iblo) = bQuickBrute(1);
        ubBlo(iblo) = bQuickBrute(2);
    end
end

elim = zeros(1,nblo);

for iblo=1:nblo
    elim(ubBlo<lbBlo(iblo))=1;
end

newblocks = blocks;
newblocks(elim==1) = [];

bounds = [lbBlo',ubBlo'];

end