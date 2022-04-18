function schemeData = solutionScheme(sadata)

%% Get Structural Analysis data
n = DAESAgetSize(sadata);
% fcn = getDAEfhandle(sadata);
sigma = DAESAgetSigma(sadata);
[c, d] = DAESAgetOffsets(sadata);
[p, q, ~, r] = DAESAgetBTF(sadata);
[~, alpha] = DAESAgetQLdata(sadata);
% constr = getConstr(sadata);
variable = DAESAgetInitData(sadata); 
alpha_eqn = getAlphaEqn(sadata);

% showStructDM(sadata);

%% Do permutations
cp = c(p); dq = d(q);
% clp = cl(p); dlq = dl(q);
sigmap = sigma(p, q);
% constrp = constr(p);
% variableq = variable(q);
alpha_eqn = alpha_eqn(p);

%% Initialize schemeData
k_d = -max(d);
schemeData(1:1-k_d) = scheme();

for k = k_d:0
    %{
    eqnp_order = cp + k;
    eqnp = find(eqnp_order >= 0);
    orderEqn = eqnp_order(eqnp);
    [block, blkposn, eqns, vars, prov, iv, solver] ...
        = analyzeStage(eqnp, orderEqn, sigmap, p, q, r, n, variable, alpha);
    %}
    
    varq_order = dq + k;
    varq = find(varq_order >= 0);
    orderVar = varq_order(varq);
    [block, blkposn, eqns, vars, prov, iv, solver, lin] ...
        = analyzeStage(varq, orderVar, cp+k, sigmap, p, q, r, n, ...
          variable, alpha, alpha_eqn);    
    
    
    i = k-k_d+1;
    
    schemeData(i) = setData(schemeData(i), k, block, blkposn, lin, eqns, ...
                                  vars, iv, prov, solver, sadata);
   
end
end

function [blku, blkposn, eqns, vars, prov, iv, solver, lin] = ...
            analyzeStage(var, orderVar, orderEqn, sigmap, p, q, r, siz, ...
            variable, alpha, alpha_eqn)
        
    onetosiz = 1:siz;
    n = length(var);
    blk = zeros(1, n);
    for i=1:n
        blk(i) = find(var(i) - r<0, 1) - 1 ; % locate the block for each variable
    end
    blku = sort(unique(blk), 'descend'); % all blocks in descending order
    
    n = length(blku); % number of blocks
    eqns = cell(n, 1);
    vars = cell(n, 1);
    prov = cell(n, 1);
    blkposn = cell(n, 1);
    iv = cell(n, 1);
    solver = zeros(n, 1);
    lin = zeros(n, 1);
    
    for i = 1:n % analyze block by block, starting from the right bottom one
        
        curblk = blku(i); % current block
        
       %% Equations to be solved
        % rc \/ diffrc = 1:siz  
        blkrg = r(curblk):r(curblk+1)-1; % block range, indices of variables in block i
        blkrest = setdiff(onetosiz, blkrg); % indices of variables outside block i
        len = r(curblk+1) - r(curblk); % block size
        blkposn{i} = blkrg; % block i's position
        
        orderEqnBlk = orderEqn(blkrg); % order of equations in block i
        inx = orderEqnBlk >= 0;
        suminx = sum(inx); % number of equations to solve in block i
        
        if suminx > 0 % if there are equations to solve
        
            eqninx = blkrg(inx); % equations to be solved in block i
            eqns{i} = [p(eqninx); orderEqn(eqninx)]'; 
            % permute back, mark down the indices of original equations and their orders

            eqnvar = max(sigmap(eqninx,:) + repmat((orderEqnBlk(inx))',[1 siz]), [], 1); % highest order of variables in these equations
            eqnvarrc = eqnvar(blkrg); % highest order of all variables in this block (not appear, -inf's)
            rc = blkrg(eqnvarrc > -1); % rc = blkrg(isfinite(eqnvarrc)), find which variables appear
            eqnvarrc = eqnvar(rc); % highest order of variables in this block that appear

            othervar = eqnvar(blkrest); % highest order of all variables outside the block (not appear, -inf's)
            diffrc = blkrest(othervar > -1); % find which variables appear
            eqnvardiffrc = eqnvar(diffrc); % mark them down      

           %% Variables to be provided
            prov{i} = [];
            m = length(rc);
            for j = 1:m
                prov{i} = [prov{i}; [repmat(q(rc(j)), [eqnvarrc(j) 1]), (0:eqnvarrc(j)-1)']];
            end
            m = length(diffrc);
            for j = 1:m
                prov{i} = [prov{i}; [repmat(q(diffrc(j)), [eqnvardiffrc(j)+1 1]), (0:eqnvardiffrc(j))']];
            end
        end
        
       %% Variables to be determined
        inx = logical((var>=r(curblk)) .* (var<r(curblk+1)));
        % should optimize here!!! use find
        
        % indices of variables in block i
        vars{i} = [q(var(inx)); orderVar(inx)]';
        % permute back, mark down the indices of original variables and their orders
        
       %% Initial values
        m = sum(inx); % number of variables in block i
        for j = 1:m
            if variable(vars{i}(j,1)) > vars{i}(j,2) 
				iv{i} = [iv{i}; vars{i}(j,:)]; % initial guess required
            end
        end
        
       %% Solver
        if suminx == 0
            solver(i) = -2; % solve nothing, simply initialize the variable
            
        elseif suminx == len % if all equation to be solved in the block       
            if alpha(blku(i)) == 1 
                solver(i) = -1;
        % Linear block, call linear solver
            elseif alpha(blku(i)) == 0 && min(orderEqnBlk)>0
                solver(i) = -1;
        % Non-linear block, minimal order of equations are larger than zero, call linear solver
            end
        end % otherwise call non-linear solver        
        
       %% Linearity
        inx = find(orderEqnBlk==0);
        if isempty(inx) % all are differentiated equations
            lin(i) = 1; % then linear
        else
            inx = blkrg(inx);
            if all(alpha_eqn(inx)==1) % alpha_eqn(i) = 1 --> f_i is linear 
                                      % in the leading derivatives only within itself
                % each equation is linear in the leading derivatives only within the equation itself
                lin(i) = 1;
            end
        end
    end
end