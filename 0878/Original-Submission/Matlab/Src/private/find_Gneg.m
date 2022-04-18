function [Gcol, Gcold] = find_Gneg(A, B, C, n, Cd)
% FIND_GNEG  Make Gj matrices for negative j
%
%   Gcol = find_Gneg(A, B, C, n) returns the column matrix
%
%          [G(-(n-p-1)); ... ; G(-2) ; G(-1)]
%
%   where G(-j) = cov(x(t),y(t+j)). [Gcol,Gcold]=find_Gneg(A,B,C,n,Cd) finds
%   also the derivative of Gcol

  DIFF = nargout > 1;
  r = size(C{1},1);
  p = size(A,2)/r;
  q = size(B,2)/r;
  m = n-p+q;
  % CHANGE A FROM [A1..Ap] to [Ap..A1]
  A = cell2mat(fliplr(makecell(A))); 
  if isempty(A), A = zeros(r,0); end
  if isempty(B), B = zeros(r,0); end
  % FIND Ccol = [C0;...; C(q); C(q+1);...; C(n-p+q)]
  Ccol = cat(1, C{:}, zeros(r*(n-p-1),r));
  Gcol = zeros(r*(n-p-1), r);
  if DIFF
    nPar = size(Cd{1},3);
    Ccold = cat(1, Cd{:}, zeros(r*(n-p-1), r, nPar));
    Gcold = zeros(r*(n-p-1), r, nPar);
  end
  for j = q+1:m-1
    J = j*r+1 : (j+1)*r;
    K = max(0,j-p)*r+1 : j*r;
    if j<p, AA = A(:,r*(p-j)+1:end); else AA = A; end
    Ccol(J,:) = AA*Ccol(K,:);
    if DIFF, Ccold(J,:,:) = ABGdiff(AA, -1, Ccol(K,:), Ccold(K,:,:)); end
  end
  % TRANSPOSE EACH C{j}:
  Ccol = reshape(permute(reshape(Ccol,r,m,r), [3,2,1]), r*m,r);
  if DIFF
    Ccold = reshape(permute(reshape(Ccold, r,m,r,nPar), [3,2,1,4]), r*m,r,nPar);
  end
  for j=1:n-p-1
    J = (n-p-j-1)*r+1 : (n-p-j)*r;
    I = j*r+1 : (j+1)*r;
    K = (j+1)*r+1 : (j+q+1)*r;
    Gcol(J,:) = Ccol(I,:) + B*Ccol(K,:);
    if DIFF
      Gcold(J,:,:) = Ccold(I,:,:) + ABGdiff(B, p*r^2+1, Ccol(K,:),Ccold(K,:,:)); 
    end
    J = J-r;
  end
end

function Fd = ABGdiff(AB, iParmat, C, Cd)
  % Differentiate A·C or B·C. iParmat = 1 for A, p·r^2+1 for B
  % Set iParmat = -1 to differentiate [Ap...A1]·C.
  nPar = size(Cd, 3);
  r = size(AB,1);
  p = size(AB,2)/r;
  Fd = reshape(AB*reshape(Cd,r*p,r*nPar),r,r,nPar);
  k = abs(iParmat);
  I = 0:p-1;
  if iParmat<0, I = fliplr(I); end
  for i = I
    for c=r*i+1:r*i+r
      for l=1:r
        Fd(l,:,k) = Fd(l,:,k) + C(c,:);
        k = k+1;
      end
    end
  end
end

