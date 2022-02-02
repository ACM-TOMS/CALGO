%MDS_ADD_PROD  Add product to matrix-function and update derivative
%
%  FSNEW = MDS_ADD_PROD(FS,X,Y,kX,kY) calculates F + X·Y and its derivatives
%  with respect to all the parameter matrices of a VARMA model, [A1...Ap B1...Bq
%  Sig]. X should be the kX-th parameter matrix and Y the kY-th one. FS is a
%  structure with the following fields:
%    mat    The r×r matrix F 
%    der    A (p+q+1)-element cell-array with the derivatives of F w.r.t. to the
%           parameter matrices; der{k} is an r×r cell matrix with derivatives
%           with respect to the k-th parameter matrix, P, and der{k}{l,c} is an
%           r×r matrix with the derivative of F with respect to P(l,c).
%    spcod  Vector of sparsity codes. The sparsity code has one of the values:
%             'z'  FS.der{k} is all zeros
%             'e'  FS.der{k}{i,j}(i,j) is the only nonzero in FS.der{k}{i,j}
%             'r'  The nonzeros of FS.der{k}{i,j} are all in its i-th row
%             'c'  The nonzeros of FS.der{k}{i,j} are all in its j-th column
%             'x'  The nonzeros of FS.der{k}{i,j} are all in its i-th column
%             'f'  All cells in FS.der{k} are full matrices
%  F + X·Y and its derivative are returned in FSNEW which is a structure as FS.
%
%  FSNEW = MDS_ADD_PROD(FS,X,Y,'T',kX,kY) calculates F + X·Y' and its
%  derivatives. FS, X and Y are as described above.
%
%  FSNEW = MDS_ADD_PROD(FS,X,A,kX) calculates F + X·A and its derivaties for a
%  constant matrix A. In this case A need not have r columns.
%
%  FSNEW = MDS_ADD_PROD(FS,X,GS,kX) calculates F + X·G and its derivatives. FS
%  and GS are structures as described above and X is the kX-th parameter matrix.
%  G need not have r columns (and if it does not it must have sparsity code 'z',
%  'r' or 'f', and now der{k} is r×size(G,2)).
%
%  FSNEW = MDS_ADD_PROD(FS,X,GS,'T',kX) calculates F + X·G' and its derivatives.
%  FS, GS and X are as above. G need not have r rows, and if it does not, must
%  have sparsity code 'f' or 'z', and now der{k} is r×size(G,1).

function FSNEW = mds_add_prod(FS,X,Y,varargin)
  N = length(FS.der);
  if ischar(varargin{1}), transp='T'; varargin(1) = []; else transp=''; end
  kX = varargin{1};
  if isstruct(Y)
    GS = Y;
    % global cover
    % for k=1:N
    %   i = 2 + (k==kX) + 2*isequal(transp,'T');
    %   j = 2 + strfind('zercxf',GS.spcod(k));
    %   cover(i,j)='x';
    % end
    cod = ['FpXG' transp];
  elseif nargin >= 5
    cod = ['FpXY' transp]; 
    kY = varargin{2};
    ascertain(kX~=kY);
  else
    cod = 'FpXA';
    kY = 0;
  end
  r = size(X,1);
  der = FS.der;
  spcod = FS.spcod;
  % fprintf('Cod=%s kX=%d\n', cod, kX);
  switch cod
    case 'FpXA'
      F = FS.mat + X*Y;
      for c = 1:r
        for l = 1:r
          der{kX}{l,c}(l,:) = der{kX}{l,c}(l,:) + Y(c,:);
        end
      end
      spcod(kX) = setsp(spcod(kX),'r');
    case 'FpXY'
      F = FS.mat + X*Y;
      for c = 1:r
        for l = 1:r
          der{kX}{l,c}(l,:) = der{kX}{l,c}(l,:) + Y(c,:);
          der{kY}{l,c}(:,c) = der{kY}{l,c}(:,c) + X(:,l);
        end
      end
      spcod(kX) = setsp(spcod(kX),'r');
      spcod(kY) = setsp(spcod(kY),'c');
    case 'FpXYT'
      F = FS.mat + X*Y';
      for c = 1:r
        for l = 1:r
          der{kX}{l,c}(l,:) = der{kX}{l,c}(l,:) + Y(:,c)';
          der{kY}{l,c}(:,l) = der{kY}{l,c}(:,l) + X(:,c);
        end
      end
      spcod(kX) = setsp(spcod(kX),'r');
      spcod(kY) = setsp(spcod(kY),'x');
    case 'FpXG'
      F = FS.mat + X*GS.mat;
      Gd = GS.der;
      for k = 1:N
        % fprintf('GS.spcod(%d)=%s\n',k,GS.spcod(k))
        for c = 1:r
          for l = 1:r
           switch GS.spcod(k)
             case 'f', der{k}{l,c} = der{k}{l,c} + X*Gd{k}{l,c};
             case 'r', der{k}{l,c} = der{k}{l,c} + X(:,l)*Gd{k}{l,c}(l,:);
             case 'c', der{k}{l,c}(:,c) = der{k}{l,c}(:,c) + X*Gd{k}{l,c}(:,c);
             case 'x', der{k}{l,c}(:,l) = der{k}{l,c}(:,l) + X*Gd{k}{l,c}(:,l);
             case 'e', der{k}{l,c}(:,c) = der{k}{l,c}(:,c) + X(:,l);
             case 'z'  % do nothing
           end
          end
        end
        switch GS.spcod(k)
          case {'f','r'}, spcod(k) = setsp(spcod(k),'f');
          case {'c','e'}, spcod(k) = setsp(spcod(k),'c');
          case {'x'},     spcod(k) = setsp(spcod(k),'x');
        end
      end
      for c = 1:r
        for l = 1:r
          der{kX}{l,c}(l,:) = der{kX}{l,c}(l,:) + GS.mat(c,:);
        end
      end
      spcod(kX) = setsp(spcod(kX),'r');
    case 'FpXGT'
      F = FS.mat + X*GS.mat';
      Gd = GS.der;
      for k = 1:N
        % fprintf('GS.spcod(%d)=%s\n',k,GS.spcod(k))
        for c = 1:r
          for l = 1:r
            switch GS.spcod(k)
              case 'f', der{k}{l,c} = der{k}{l,c} + X*Gd{k}{l,c}';
              case 'c', der{k}{l,c} = der{k}{l,c} + X(:,c)*Gd{k}{l,c}(:,c)';
              case 'x', der{k}{l,c} = der{k}{l,c} + X(:,l)*Gd{k}{l,c}(:,l)';
              case 'r', der{k}{l,c}(:,l) = der{k}{l,c}(:,l)+X*Gd{k}{l,c}(l,:)';
              case 'e', der{k}{l,c}(:,l) = der{k}{l,c}(:,l)+X(:,c);
              case 'z'  % do nothing
            end
          end
        end
        switch GS.spcod(k)
          case {'r','e'},     spcod(k) = setsp(spcod(k),'x');
          case {'f','c','x'}, spcod(k) = setsp(spcod(k),'f');
          case {'z'},         % do nothing
        end
      end
      for c = 1:r
        for l = 1:r
          der{kX}{l,c}(l,:) = der{kX}{l,c}(l,:) + GS.mat(:,c)';
        end
      end
      spcod(kX) = setsp(spcod(kX),'r');
  end
  FSNEW.mat = F;
  FSNEW.der = der;
  FSNEW.spcod = spcod;
end

function spcod = setsp(spcod, s)
  % If s~=spcod and either both are in rcx or both are in ex then spcod:='f'. 
  % Otherwise if s is >spcod in the ordering z<e<r<c<f then spcod:=s
  if spcod ~= s && ...
    (any('xcr'==spcod) && any('xcr'==s) || any('xe'==spcod) && any('xe'==s))
      spcod = 'f';
  elseif strfind('zexcrf',s) > strfind('zexcrf',spcod)
    spcod = s;
  end
end
