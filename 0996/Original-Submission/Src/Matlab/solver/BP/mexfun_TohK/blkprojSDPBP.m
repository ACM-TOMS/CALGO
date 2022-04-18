%%***********************************************************************
%% compute projection onto the cone of 
%% positive semidefinite matrices
%%
%%***********************************************************************
%% SDPNAL+ 
%% Copyright (c) 2014 by
%% Liuqin Yang, Defeng Sun, and Kim-Chuan Toh
%%***********************************************************************

    function [Xp,d,rranknew] = blkprojSDPBP(X,K1,rrank)
    
    if (nargin < 3); rrank = []; end
    optseig.jobz = 1; optseig.range = 2; 
    optseig.abstol = 1e-16; optseig.eigchoice = 1;
    can_use_eigs = 1;     
    smtol = 0; 
    ranktol = 1e-15; 
%%    
    [rowSize,colSize] = size(X);     
    if isfield(K1,'l') && (K1.l > 0)
       nDim = K1.l + sum(K1.s .* K1.s); 
       mDim = K1.l + sum(K1.s); 
       Xp = zeros(rowSize,colSize);
       d = zeros(mDim,1); 
       Xp(1:K1.l) = max(X(1:K1.l),0);
       d(1:K1.l) = X(1:K1.l); 
       colPointer = K1.l; 
       eigenValpointer = K1.l; 
    else
       nDim = sum(K1.s .* K1.s); 
       mDim = sum(K1.s); 
       Xp = zeros(rowSize,colSize);
       d = zeros(mDim,1); 
       colPointer = 0;        
       eigenValpointer = 0; 
    end
    rranknew = zeros(1,length(K1.s));
%%    
    if (length(K1.s)==1) %%single SDP block
       ns = K1.s;         
       if (colPointer==0)
          XMat = reshape(X,ns,ns);
       else
          XMat = reshape(X(colPointer+[1:K1.s*K1.s]),ns,ns); 
       end
       use_eigs = 0;
       flag_eigs = 0; 
       if isempty(rrank)
          rrank = ns/2; can_use_eigs = 0; 
       end       
       if (can_use_eigs) && (ns >= 2000) && (rrank < min(0.1*ns,20))
          use_eigs = 1; 
       elseif (can_use_eigs) && (ns >= 2000) && (ns-rrank < min(0.1*ns,20))
          use_eigs = 2; 
       end 
       if (use_eigs==1) 
          [psdXMat,eigXMat,flag_eigs] = projSDPsub_eigs(XMat,rrank+5,smtol); 
          rranknew = length(find(eigXMat > ranktol*max(abs(eigXMat)))); 
          fprintf('+')
       elseif (use_eigs==2) 
          [Xtmp,dtmp,flag_eigs] = projSDPsub_eigs(-XMat,ns-rrank+5,smtol); 
          psdXMat = XMat+Xtmp; 
          rranknew = ns-length(find(dtmp > ranktol*max(abs(dtmp)))); 
          fprintf('#')
       end
       if (flag_eigs) %% recompute if eigs fails.
          use_eigs = 0;  
          fprintf('-')
       end 
       if (use_eigs==0)
          if (rrank <= ns/2)
             sign = 1;
          else
             sign = -1;
          end 
          if (optseig.eigchoice==1) ...
             && ((rrank >= 0.85*ns) || (rrank > 0 && rrank <= 0.05*ns))
             [psdXMat,eigXMat] = projSDPsubpartial(XMat,smtol,sign,optseig);
             if (sign == 1)
                rranknew = length(find(eigXMat > ranktol*max(abs(eigXMat)))); 
             else                                                           
                rranknew = ns - length(find(eigXMat > ranktol*max(abs(eigXMat))));
             end
          else
             [psdXMat,eigXMat] = projSDPsub(XMat,smtol);
             rranknew = length(find(eigXMat > ranktol*max(abs(eigXMat))));
          end                   
       end
       if (colPointer==0)
          Xp = reshape(psdXMat,1,ns*ns);
       else
          Xp(colPointer+[1:ns*ns])  = reshape(psdXMat,1,ns*ns); 
       end
       if use_eigs == 2
           d(eigenValpointer+[1:ns]) = 0; % need to be fixed
       else
           d(eigenValpointer+[1:ns]) = eigXMat; 
       end
    elseif (length(K1.s)>1) %% multiple SDP blocks    
        for k=1:length(K1.s)
           nk = K1.s(k); 
           XMat = reshape(X(colPointer+[1:nk*nk]),nk,nk);
           if (nk > 1)
              [psdXMat,eigXMat] = projSDPsub(XMat,smtol);
              Xp(colPointer+[1:nk*nk])  = reshape(psdXMat,1,nk*nk); 
              d(eigenValpointer+[1:nk]) = eigXMat; 
           elseif (nk == 1) 
              Xp(colPointer+1) = max(XMat,0); 
              d(eigenValpointer+1) = XMat; 
           end
           colPointer = colPointer + nk*nk; 
           eigenValpointer = eigenValpointer + nk; 
        end
    end                                        
%%********************************************************************
%% compute projection on PSD cone using the mexeig function. 
%%********************************************************************
    function [Xp,eigX,V] = projSDPsub(X,smtol)
        
    n = length(X);
    [V,D] = mexeigK(full(X));
    eigX = diag(D);
    posidx = find(eigX > smtol);
    len = length(posidx); 
    if isempty(posidx)
       Xp = sparse(n,n); 
    elseif (length(posidx) == n)
       Xp = X; 
    else
       if (length(posidx) <= n/2) 
          dp = abs(eigX(posidx)); 
          Vtmp = V(:,posidx)*sparse(1:len,1:len,sqrt(dp));%spdiags(sqrt(dp),0,len,len); 
          Xp = Vtmp*Vtmp'; Xp = 0.5*(Xp+Xp'); 
       else
          negidx = find(eigX <= smtol);
          lenneg = length(negidx); 
          dn = abs(eigX(negidx)); 
          Vtmp = V(:,negidx)*sparse(1:lenneg,1:lenneg,sqrt(dn));%spdiags(sqrt(dn),0,lenneg,lenneg); 
          Xn = Vtmp*Vtmp'; Xn = 0.5*(Xn+Xn'); 
          Xp = X+Xn;
       end
    end
%%********************************************************************
%% compute projection on PSD cone using the mexDSYEVX function. 
%%********************************************************************
   function [Xp,eigX,V] = projSDPsubpartial(X,smtol,sign,optseig)
        
    if (nargin < 4)
       optseig.jobz= 1; optseig.range = 2; optseig.abstol = 1e-16; 
    else
       if ~isfield(optseig,'jobz'); optseig.jobz= 1; end
       if ~isfield(optseig,'range'); optseig.range = 2; end
       if ~isfield(optseig,'abstol'); optseig.abstol = 1e-16; end
    end    
    if (sign == 1)
       [V,D] = mexeigPartialK(full(X), smtol, 1e15, optseig); 
       eigX = D;
       nD = (D>smtol);   
       Vtmp = V(:,nD);
       Xp = (Vtmp*diag(D(nD))*Vtmp');
    else
       [V,D] = mexeigPartialK(full(-X), smtol, 1e15, optseig);
       eigX = D;
       nD = (D>smtol);   
       Vtmp = V(:,nD);
       Xn = (Vtmp*diag(D(nD))*Vtmp');
       Xp = (Xn+X);
    end
    Xp = 0.5*(Xp+Xp');
%%********************************************************************
%% compute projection onto PSD cone using the eigs function. 
%%********************************************************************
    function [Xp,eigX,V,flag] = projSDPsub_eigs(X,rrank,smtol)

    n = length(X);
    options.disp  = 0; 
    %options.issym = 1; % not necessary since 'LA' is specified 
    options.tol   = 1e-14; 
    [V,D,flag] = eigs(full(X),rrank,'LA',options);
    eigX = diag(D);
    posidx = find(eigX > smtol);
    len = length(posidx); 
    if isempty(posidx)
       Xp = sparse(n,n); 
    elseif (length(posidx) == n)
       Xp = X; 
    else
       dp = abs(eigX(posidx)); 
       Vtmp = V(:,posidx)*spdiags(sqrt(dp),0,len,len); 
       Xp = Vtmp*Vtmp'; Xp = 0.5*(Xp+Xp'); 
    end
%%********************************************************************
