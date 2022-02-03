function [IVset, constraintSet, alphaC, alphaF, linF] = ...
    QuasiLinAnal(fcn, n, sigmaPerm, c, ccl, cfl, dfl, ...
    p, q, rC, rF, JNZ, varargin)


%% Initialization
numFineBlk = length(rF)-1; % number of fine blocks
numCoarseBlk = length(rC) -1; % number of coarse blocks

global DAEZERO;

% Permute offsets
c = c(p);
ccl = ccl(p);
cfl = cfl(p);
dfl = dfl(q);
JNZ = JNZ(p,q);

% Initialize output
IVset = zeros(1,n);
constraintSet = zeros(n,1);
alphaF = true(1,numFineBlk);
alphaC  = true(1,numCoarseBlk);

%% QLA for fine BTF
t = qla(n); t = setType(t, 1:n, zeros(1, n));

DAEZERO = qla(n);
y(1:n) = t; y=y';

for i = 1:numFineBlk
    rc = rF(i):rF(i+1)-1; % block range
    qrc = q(rc);
    prc = p(rc);
    len = rF(i+1) - rF(i); % block size
    sigmaBlock = sigmaPerm(rc, rc);
    
    % Only checks the entries with equality d_j-c_i=s_ij
    JNZBlock = JNZ(rc, rc);
    ind = JNZBlock==0;
    sigmaBlock(ind) = inf;
    
    for j = 1:len
        index = qrc(j); % index of variable
        sigmaj = sigmaBlock(:,j);
        
        y(index) = setOffsets(y(index), prc, sigmaj);
        % set offsets for the equations with indices prc
        
        type = zeros(1,len);
        indx = sigmaj==0;
        type(indx) = 2; % linear
        indx = sigmaj>0;
        type(indx) = 1; % immediate
        
        y(index) = setType(y(index), prc, type);
    end
end

f = feval(fcn, t, y, varargin{:}); % Evaluate the DAE
indx = [];
for i=1:n
    if getType(f(i),i)==3
        indx = [indx i]; % mark the NQL eqns
    end
end
linF=~sparse(indx,1,1,n,1); % assign ones to the QL eqns

for i=1:numFineBlk
    rc = rF(i):rF(i+1)-1;
    for j = rc
        if linF(p(j))==0 && cfl(j)==0
            % use local fine offsets to check if ith fine block is NQL
            alphaF(i) = false;
            break;
        end
    end
    IVset(rc) = dfl(rc) - alphaF(i) + 1;
    constraintSet(rc) = c(rc) - alphaF(i) + 1;
end

% Permute back
constraintSet(p) = constraintSet;
IVset(q) = IVset;

%% QLA for coarse BTF
t = qla(n); t = setType(t,1:n, zeros(1, n));

DAEZERO = qla(n);
y(1:n) = t; y=y';

for i = 1:numCoarseBlk
    rc = rC(i):rC(i+1)-1; % block range
    qrc = q(rc);
    prc = p(rc);
    len = rC(i+1) - rC(i); % block size
    sigmaBlock = sigmaPerm(rc, rc);
    
    % Only checks the entries with equality d_j-c_i=s_ij
    JNZBlock = JNZ(rc, rc);
    ind = JNZBlock==0;
    sigmaBlock(ind) = inf;
    
    for j = 1:len
        index = qrc(j); % index of variable
        sigmaj = sigmaBlock(:,j);
        
        y(index) = setOffsets(y(index), prc, sigmaj);
        % set offsets for the equations with indices prc
        
        type = zeros(1,len);
        
        indx = sigmaj==0;
        type(indx) = 2; % linear
        indx = sigmaj>0;
        type(indx) = 1; % immediate
        
        y(index) = setType(y(index), prc, type);
    end
end

f = feval(fcn, t, y, varargin{:}); % Evaluate the DAE

for i=1:numCoarseBlk
    rc = rC(i):rC(i+1)-1;
    for j = rc
        if getType(f(p(j)),p(j))==3 && ccl(j)==0
            % use local fine offsets to check if ith fine block is NQL
            alphaC(i) = false;
            break;
        end
    end
end
end