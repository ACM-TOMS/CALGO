function [x,y,info] = mosekSedumi(A,b,c,K)
%MOSEKSEDUMI  Solve the cone problem in sedumi format by MOSEK
%
%    CAUTION: Currently this function does not support Lorentz cone.
%
%    Usage:
%       [x,y,info] = mosekSedumi(A,b,c,K);

    if isfield(K,'q')
        error('MOSEKSEDUMI does not support Lorentz cone.');
    end
    if size(c,1)~=1; c = c'; end
    if size(b,2)~=1; b = b'; end
    A = A';
    %mm = size(A,1);
    
    if ~isfield(K,'f'); K.f = 0; end
    if ~isfield(K,'l'); K.l = 0; end
        

    nn = K.f + K.l;
    prob.a = A(:,1:nn); AA = A(:,nn+1:end);%;sparse(mm,nn);
    prob.c = c(1:nn)'; cc = c(nn+1:end);%sparse(1,nn);
    prob.blx = [-inf*ones(K.f,1); zeros(K.l,1)];
    
    
    pp = length(K.s);
    dim2 = K.s .* K.s;
    pointer = cumsum([0,dim2]);
    

    
    prob.bardim = K.s;

    barc.subj = []; barc.subk =[]; barc.subl = []; barc.val = [];
    for jj = 1:pp
        [~,row,col,val] = vecToTrilMat( cc( (1:dim2(jj)) + pointer(jj) ) );
        barc.subj = [barc.subj; jj*ones(length(row),1)];
        barc.subk = [barc.subk; row];
        barc.subl = [barc.subl; col];
        barc.val  = [barc.val;  val];
    end
    prob.barc = barc;
    
    bara.subi = []; bara.subj = []; bara.subk =[]; bara.subl = []; bara.val = [];
    for jj = 1:pp
        [ii,kk,ll,val] = vecToTrilMat( AA( :, (1:dim2(jj)) + pointer(jj) ) );
        bara.subi = [bara.subi; ii];
        bara.subj = [bara.subj; jj*ones(length(ii),1)];
        bara.subk = [bara.subk; kk];
        bara.subl = [bara.subl; ll];
        bara.val  = [bara.val; val];
    end
    prob.bara = bara;

    prob.blc = b;
    prob.buc = b;

    param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-8;
    [~,res] = mosekopt('minimize info', prob, param);
    
    eidx = cumsum(prob.bardim.*(prob.bardim+1)/2);
    sidx = [0,eidx(1:end-1)]+1;
    x = res.sol.itr.xx;
    for ii = 1:pp
        x = [x; trilVecToVecMat(res.sol.itr.barx(sidx(ii):eidx(ii)))];
    end
    y = res.sol.itr.y;
    info.cputime = res.info.MSK_DINF_OPTIMIZER_TIME;
end



function [ii,kk,ll,val] = vecToTrilMat(vec)
    dim = sqrt(size(vec,2));
    [ii, idx, val] = find(vec);
    ll = rem( idx-1, dim ) + 1;
    kk = floor( (idx-1) / dim ) + 1;
    trilidx = kk >= ll;
    ii  = ii(trilidx);
    kk  = kk(trilidx);
    ll  = ll(trilidx);
    val = val(trilidx);
end

function vecMat = trilVecToVecMat(trilVec)
    dim = sqrt(length(trilVec)*2+0.25)-0.5;
    [row,col] = find(tril(ones(dim)));
    mat = sparse(row,col,trilVec);
    mat = mat + tril(mat,-1)';
    vecMat = mat(:);
end