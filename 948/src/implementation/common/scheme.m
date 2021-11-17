classdef scheme
    
    properties(GetAccess = 'private', SetAccess = 'private')
        k; blk; blkposn; lin; eqn; var; iv; prov; solver; sadata;
    end
    
    methods(Hidden=true)
        % methods, including the constructor are defined in this block
        function schemeData = scheme(k, blk, blkposn, lin, eqn, var, iv, prov, solver, sadata)
            
            % schemeData contains max(d)+1 'scheme' objects. Each includes
            %
            % k:
            %    an integer for solving stage
            %    n refers the number of blocks that should be solved in stage k
            %
            % blk:
            %    array, length of n, blocks to be solved in this stage
            %
            % blokposn:
            %    array, block positions
            %
            % lin:
            %    array, length of n, quasi-linearity of each block
            %
            % eqn:
            %    cell, length of n, each cell gives which equations to be solved in the block
            %
            % var:
            %    cell, length of n, each cell gives which variables should be
            %    determined in the block
            %
            % iv:
            %    cell, length of n, each cell indicates which variables need to be initialized
            %    in this block
            %
            % prov (should be renamed):
            %    cell, length of n, each cell indicates which variables already solved should
            %    be provided in order to solve this block
            %
            % solver:
            %    array, length of n, call linear solver if -1; call nonlinear solver if 0;
            
            if nargin == 0
                schemeData.k = NaN;
                schemeData.blk = NaN;
                schemeData.blkposn = NaN;
                schemeData.lin = NaN;
                schemeData.eqn = NaN;
                schemeData.var = NaN;
                schemeData.iv = NaN;
                schemeData.prov = NaN;
                schemeData.solver = NaN;
                schemeData.sadata = NaN;
            else
                schemeData.k = k;
                schemeData.blk = blk;
                schemeData.blkposn = blkposn;
                schemeData.lin = lin;
                schemeData.eqn = eqn;
                schemeData.var = var;
                schemeData.iv = iv;
                schemeData.prov = prov;
                schemeData.solver = solver;
                schemeData.sadata = sadata;
            end
        end
        
        function schemeData = setData(schemeData, k, blk, blkposn, lin, ...
                eqn, var, iv, prov, solver, sadata)
            schemeData.k = k;
            schemeData.blk = blk;
            schemeData.blkposn = blkposn;
            schemeData.lin = lin;
            schemeData.eqn = eqn;
            schemeData.var = var;
            schemeData.iv = iv;
            schemeData.prov = prov;
            schemeData.solver = solver;
            schemeData.sadata = sadata;
        end
        
        function [k, blk, blkposn, lin, iv, prov, eqn, var, solver, ...
                sadata] = getAll(schemeData)
            k = schemeData.k;
            blk = schemeData.blk;
            blkposn = schemeData.blkposn;
            lin = schemeData.lin;
            iv = schemeData.iv;
            prov = schemeData.prov;
            eqn = schemeData.eqn;
            var = schemeData.var;
            solver = schemeData.solver;
            sadata = schemeData.sadata;
        end
        
        function [sadata] = getSAdata(schemeData)
            sadata = schemeData.sadata;
        end
    end
end