classdef SAdata

    properties(GetAccess = 'private', SetAccess = 'private', Hidden = true)
        fcn;
        n;
        sigma;
        DOF;
        HVT;
        index;
        c; d;     % global       offsets
        ccl; dcl; % coarse local offsets
        cfl; dfl; % fine   local offsets   
        IVset; constraintSet;
        p; q; r; rf;
        alpha; alphaf;
        alpha_eqn;
        JNZ;
        SWP;
        exitflag;
        time;
        missEqns; missVars;
        param;
    end
    
    methods(Hidden = true)
        % methods, including the constructor are defined in this block
        function sadata = SAdata(n, fcn, sigma, DOF, HVT, index, c, d, ...
                ccl, dcl, cfl, dfl, IVset, constraintSet, p, q, r, rf, ...
                alpha, alphaf, alpha_eqn, JNZ, SWP, exitflag, time, ...
                missEqns, missVars, param)
            
            if nargin == 1
                sadata.fcn = blanks(0);
                sadata.n = NaN;
                sadata.sigma(1:n, 1:n) = NaN;
                sadata.DOF = NaN;
                sadata.HVT(1:n) = NaN;
                sadata.index = NaN;
                sadata.c(1:n) = NaN;
                sadata.d(1:n) = NaN;
                sadata.ccl(1:n) = NaN;
                sadata.dcl(1:n) = NaN;
                sadata.cfl(1:n) = NaN;
                sadata.dfl(1:n) = NaN;
                sadata.IVset= NaN;
                sadata.constraintSet = NaN;
                sadata.p = NaN;
                sadata.q = NaN;
                sadata.r = NaN;
                sadata.rf = NaN;
                sadata.alpha = NaN;
                sadata.alphaf = NaN;
                sadata.alpha_eqn = NaN;
                sadata.JNZ = NaN;
                sadata.SWP = NaN;
                sadata.exitflag = NaN;
                sadata.time = NaN;
                sadata.missEqns = NaN;
                sadata.missVars = NaN;
                sadata.param = NaN;
                
            else
                sadata.fcn = fcn;
                sadata.n = n;
                sadata.sigma = sigma;
                sadata.DOF = DOF;
                sadata.HVT(1:n) = HVT;
                sadata.index = index;
                sadata.c = c';
                sadata.d = d;
                sadata.ccl = ccl';
                sadata.dcl = dcl;
                sadata.cfl = cfl';
                sadata.dfl = dfl;
                sadata.IVset= IVset;
                sadata.constraintSet = constraintSet;
                sadata.p = p;
                sadata.q = q;
                sadata.r = r;
                sadata.rf = rf;
                sadata.alpha = alpha;
                sadata.alphaf = alphaf;
                sadata.alpha_eqn = alpha_eqn;
                sadata.JNZ = JNZ;
                sadata.SWP = SWP;
                sadata.exitflag = exitflag;
                sadata.time = time;
                sadata.missEqns = missEqns;
                sadata.missVars = missVars;
                
                if nargin == 28
                    sadata.param = param;
                elseif nargin == 27
                    sadata.param = [];
                end
            end
        end
        
        function n = DAESAgetSize(sadata)
            n = sadata.n;
        end
        
        function fcn = DAESAgetDAEfhandle(sadata)
            fcn = sadata.fcn;
        end
        
        function sigma = DAESAgetSigma(sadata)
            sigma = sadata.sigma;
        end
        
        function DOF = DAESAgetDOF(sadata)
            DOF = sadata.DOF ;
        end
        
        function index = DAESAgetIndex(sadata)
            index = sadata.index;
        end
        
        function [c, d, C, D] = DAESAgetOffsets(sadata)
            c = sadata.c;
            d = sadata.d;
            C = sadata.cfl;
            D = sadata.dfl;
        end
        
        function iv = DAESAgetInitData(sadata)
            iv = sadata.IVset;
        end
        
        function constr = DAESAgetConstr(sadata)
            constr = sadata.constraintSet;
        end
        
        function [ql, qlf] = DAESAgetQLdata(sadata)
            ql = sadata.alpha;
            qlf = sadata.alphaf;
        end
        
        function swp = DAESAisSWP(sadata)
            swp = sadata.SWP;
        end
        
        function [missEqns, missVars] = DAESAgetMissing(sadata)
            missEqns = sadata.missEqns;
            missVars = sadata.missVars;
        end
        
        function [pe,pv,cb,fb] = DAESAgetBTF(sadata)
            pe = sadata.p; pv = sadata.q;
            switch sadata.exitflag
                case 0
                    pe = sadata.p; pv = sadata.q;
                    fb = sadata.rf; cb = sadata.r;
                case {-1, -2};
                    cb = []; fb = [sadata.r; sadata.rf];
                case -3
                    cb = []; fb = [];
            end
        end        
        
        function HVT = DAESAgetHVT(sadata)
            HVT = sadata.HVT;
        end
        
        function disp(sadata)
            fprintf('    an SAdata object for ''%s'' function\n\n', ...
                func2str(sadata.fcn));
        end
    end
    
    methods(Hidden = true)

        
        function lin = getAlphaEqn(sadata)
            lin = sadata.alpha_eqn;
        end
        
        function JNZ = getJNZ(sadata)
            JNZ = sadata.JNZ;
        end
        
        function exitflag = getExitflag(sadata)
            exitflag = sadata.exitflag;
        end
        
        function param = getParam(sadata)
            param = sadata.param;
        end
        
        function time = getCPUtime(sadata)
            time = sadata.time;
        end
        
        function [n, fcn, sigma, DOF, HVT, index, c, d, ...
                clocal, dlocal, iv, constr, p, q, r, rf, alpha, alphaf, ...
                alpha_eqn, JNZ, swp, exitflag, time, missEqns, missVars, ...
                param] = ...
                getAll(sadata)
            
            fcn = sadata.fcn;
            n = sadata.n;
            sigma = sadata.sigma;
            DOF = sadata.DOF;
            HVT = sadata.HVT;
            index = sadata.index;
            c(1:n) = sadata.c;
            d(1:n) = sadata.d;
            clocal(1:n) = sadata.cfl;
            dlocal(1:n) = sadata.dfl;
            iv = sadata.IVset;
            constr = sadata.constraintSet;
            p = sadata.p;
            q = sadata.q;
            r = sadata.r;
            rf = sadata.rf;
            alpha = sadata.alpha;
            alphaf = sadata.alphaf;
            alpha_eqn = sadata.alpha_eqn;
            JNZ = sadata.JNZ;
            swp = sadata.SWP;
            exitflag = sadata.exitflag;
            time = sadata.time;
            missEqns = sadata.missEqns;
            missVars = sadata.missVars;
            param = sadata.param;
        end
    end
end
