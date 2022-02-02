function [x,y,sdpt3Info] = sdpt3Sedumi(At,b,c,K,OPTIONS)
%SDPT3SEDUMI   Solve the cone problem in sedumi format by SDPT3
%
%   Usage:
%      [x,y,info] = SDPT3SEDUMI(At,b,c,K);
%      [x,y,info] = SDPT3SEDUMI(At,b,c,K,OPTIONS);
    if nargin == 4
        OPTIONS.printlevel = 2;
        OPTIONS.gaptol     = 1e-8;
        OPTIONS.inftol     = 1e-8;
        OPTIONS.steptol    = 1e-6;
        OPTIONS.maxit      = 100;
    end
    if size(At,1) < size(At,2)
       At = At';  
    end    
    smallbkldim = 40;
    [blkSdpt3,AtSdpt3,CSdpt3,bSdpt3,perm] = read_sedumi(At',b,c,K,smallbkldim);
    if OPTIONS.printlevel <= 1
        fprintf('\nSDPT3: Infeasible path-following algorithms\n');
    end
    [~,XSdpt3,y,~,sdpt3Info] = sqlp(blkSdpt3,AtSdpt3,CSdpt3,bSdpt3,OPTIONS);
    %     pdipmInfo.cpusec = pdipmInfo.cputime;
    %     pdipmInfo.numerr = pdipmInfo.termcode;
    clear  blkSdpt3 AtSdpt3 CSdpt3 bSdpt3
    sdpPosVect = 0;
    sdppointer = 0;
    xLength = 0;
    if isfield(K,'f') && ~isempty(K.f) && K.f > 0
        xLength  = xLength + K.f;
    end
    if isfield(K,'l') && ~isempty(K.l) && K.l > 0
        xLength  = xLength + K.l;
    end
    if isfield(K,'s') && ~isempty(K.s)
        for i=1:length(K.s)
            xLength  = xLength + K.s(i)*K.s(i);
            sdppointer = sdppointer+K.s(i)*K.s(i);
            sdpPosVect = [sdpPosVect, sdppointer];
        end
    end
    x = zeros(xLength,1);
    noBlockSDPT3 = size(XSdpt3,1);
    nonSDPpointer = 0;
    for i=1:noBlockSDPT3
        sizeSDPT3 = size(XSdpt3{i},2);
        if sizeSDPT3 == 1
            x(nonSDPpointer+1:nonSDPpointer+size(XSdpt3{i},1),1) = XSdpt3{i};
            nonSDPpointer = nonSDPpointer+size(XSdpt3{i},1);
        else
            lenPerm = length(perm{i});
            if lenPerm == 1
                p = perm{i};
                x(nonSDPpointer+sdpPosVect(p)+1:nonSDPpointer+sdpPosVect(p+1),1) = reshape(XSdpt3{i},sizeSDPT3*sizeSDPT3,1);
            else
                blockPointer = 0;
                for j=1:lenPerm
                    p = perm{i}(j);
                    oneBlock = XSdpt3{i}(blockPointer+1:blockPointer+K.s(p),blockPointer+1:blockPointer+K.s(p));
                    x(nonSDPpointer+sdpPosVect(p)+1:nonSDPpointer+sdpPosVect(p+1),1) = ...
                        reshape(oneBlock,K.s(p)*K.s(p),1);
                    blockPointer = blockPointer + K.s(p);
                end
            end
        end
    end
    clear XSdpt3    
end

