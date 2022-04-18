function [c, d, ccl, dcl, HVT, DOF, p, q, r, cc, rr, exitflag] = lapdm(A)

[A, p, q, r, s, cc, rr] = permuteMatrix(A);
% Permute A into UTF

if any(r~=s) % cannot solve a LAP    
    % Get the coarse blocks, 4-by-4
    c=[]; d=[]; ccl=[]; dcl=[]; HVT=NaN; DOF=NaN;
    r=[]; exitflag=-1;
    return;
else
    exitflag = 0;
end

% Find HVT in each block
HVT = findHVT(A, r);

% Find global and coarse local offsets in each block using HVT
[c, d, ccl, dcl, HVT] = compOffsets(A, HVT, r);

% Degrees of freedom
DOF = sum(d) - sum(c);

% Permute global offsets back
c(p) = c;      d(q) = d;
ccl(p) = ccl;  dcl(q) = dcl;

n = length(c);
HVTmat = sparse(HVT, 1:n, ones(1,n), n, n);
HVTmat(p,q) = HVTmat;
[HVT, ~] = find(HVTmat);
end
