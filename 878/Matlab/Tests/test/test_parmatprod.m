%TEST_PARMATPROD  Test parameter matrix product derivative functions
%
%  TEST_PARMATPROD tests mds_set_zero, mds_set_parmat and mds_add_prod. Use
%  TEST_PARMATPROD QUIET to reduce output.
%  
%  For mds_add_prod(FS,X,GS,kX) and mds_add_prod(FS,X,GS,'T',kX) the tests are
%  designed to cover all possible combinations of spcod and transpose on G both
%  for derivative w.r.t. X and w.r.t. another parameter. To ascertain this, 
%  uncomment the block referring to cover both in mds_add_prod and here.
%
function diffout = test_parmatprod(verbosity)
  %   global cover
  %   cover = [
  %     '  zercxf'
  %     'A       '
  %     'X       '
  %     'AT      '
  %     'XT      '];
  verbose = (nargin < 1);
  fprintf('TESTING PARAMETER MATRIX PRODUCT DERIVATIVES...')
  fprintf_if(verbose,'\n')
  rand('state',1);
  r = 3;
  np = 3;
  n = np*r^2;
  x = 10*rand(n,1);
  A = rand(r,r);
  fprintf_if(verbose,'  Max numerical/actual relative derivative difference:\n')
  testcases={'A1','A2','A3','B1','B2','C1','C2','C3','D1','D2','D3','G1',...
    'G2','H1','H2'};
  for i = 1:length(testcases)
    d(i) = diff_test(@fun, x, r, np, A, testcases{i});
    fprintf_if(verbose, '  Testcase %s: %9.1e\n', testcases{i}, d(i));
    ascertain(d(i)<1e-8);
  end
  if nargout > 0, diffout = d; end  
  % cover
  disp('  OK');
  
end

function [f,g] = fun(x, r, np, A, tstcase)
  % returns a column m-vector f and an m×n Jacobian where n = length(x)
  x = x(:);
  P = reshape(x,r,[]);
  X = P(:,1:r);
  Y = P(:,r+1:2*r);
  Z = P(:,2*r+1:3*r);
  a = A(:,1);
  if nargout == 1
    switch tstcase
      case 'A1', F = X*Y;
      case 'A2', F = X*Y';
      case 'A3', F = X*A;
      case 'B1', F = X + Y*Z;
      case 'B2', F = X*Y + Y*Z;
      case 'C1', F = X*(Y);
      case 'C2', F = X*(Y)';
      case 'C3', F = X*(X) + X*(X)';
      case 'D1', F = X*(Y*Z') + Z*(Y*Z');
      case 'D2', F = X*(Y*X) + X*(Y*X)';
      case 'D3', F = X*(Y*Z')' + Z*(Y*Z')';
      case 'G1', F = Z + X*(Y + Y*Z + X*(X*Y+X*Y*Z))';
      case 'G2', F = Z*(X*Y + Y*Z + X*(X*Y + Y*Z)');
      case 'H1', F = X*a;
      case 'H2', F = Y*(X*[A -A]);
      otherwise, error(['unknown case ' tstcase]);
    end
    f = F(:)';
  else
    switch tstcase
      case 'A1'                         % F = X*Y
        F = setprod(np,X,Y,1,2);
      case 'A2'                         % F = X*Y'
        F = setprod(np,X,Y,'T',1,2);
      case 'A3'                         % F = X*A
        F = setprod(np,X,A,1);
      case 'B1'                         % F = X + Y*Z
        G = mds_set_parmat(np,X,1);
        F = mds_add_prod(G,Y,Z,2,3);
      case 'B2'                         % F = X*Y + Y*Z
        G = setprod(np,X,Y,1,2);
        F = mds_add_prod(G,Y,Z,2,3);
      case 'C1'                         % F = X*(Y)
        G = mds_set_parmat(np,Y,2);
        F = setprod(np,X,G,1);
      case 'C2'                         % F = X*(Y)'
        G = mds_set_parmat(np,Y,2);
        F = setprod(np,X,G,'T',1);
      case 'C3'                         % F = X*(X) + X*(X)'
        G = mds_set_parmat(np,X,1);
        F = setprod(np,X,G,1);
        F = mds_add_prod(F,X,G,'T',1);
      case 'D1'                         % F = X*(Y*Z') + Z*(Y*Z')
        G = setprod(np,Y,Z,'T',2,3);
        F = setprod(np,X,G,1);
        F = mds_add_prod(F,Z,G,3);
      case 'D2'                         % F = X*(Y*X) + X*(Y*X)'
        G = setprod(np,Y,X,2,1);
        F = setprod(np,X,G,1);
        F = mds_add_prod(F,X,G,'T',1);
      case 'D3'                         % F = X*(Y*Z')' + Z*(Y*Z')'
        G = setprod(np,Y,Z,'T',2,3);
        F = setprod(np,X,G,'T',1);
        F = mds_add_prod(F,Z,G,'T',3);
      case 'G1'                         % F = Z+Y*(Y+Y*Z+X*(X*Y+X*Y*Z))'
        G0 = setprod(np,Y,Z,2,3);
        G1 = setprod(np,X,Y,1,2);
        G1 = mds_add_prod(G1,X,G0,1);
        G2 = mds_set_parmat(np,Y,2);
        G2 = mds_add_prod(G2,Y,Z,2,3);
        G2 = mds_add_prod(G2,X,G1,1);
        F = mds_set_parmat(np,Z,3);
        F = mds_add_prod(F,X,G2,'T',1);
      case 'G2'                         % F = Z*(X*Y + Y*Z + X*(X*Y + Y*Z)');
        G = setprod(np,X,Y,1,2);
        G = mds_add_prod(G,Y,Z,2,3);
        G = mds_add_prod(G,X,G,'T',1);
        F = setprod(np,Z,G,3);
      case 'H1'                         % F = X*a;
        F = setprod(np,X,a,1);
      case 'H2'
        G = setprod(np,X,[A -A],1);
        F = setprod(np,Y,G,2);
      otherwise
        error(['unknown case ' tstcase]);
    end
    f = F.mat(:);
    nf = length(f);
    g=[];
    for k=1:np
      g = [g reshape(cat(2,F.der{k}{:}),nf,[])];
    end
  end
end

%SETPROD  Set matrix and its derivative to a product of matrices
%  FS = SETPROD(nparmat,X,Y,kX,kY) sets FS to a structure appropriate for
%  mds_add_prod with the product X·Y and its derivative w.r.t. all the nparmat
%  parameter matrices (nparmat should be p+q+1). The following can also be used,
%  with meanings as in mds_add_prod:
%    FS = SETPROD(nparmat,X,Y,'T',kX,kY) FS = SETPROD(nparmat,X,A,'T',kX)
%    FS = SETPROD(nparmat,X,G,kX)
%    FS = SETPROD(nparmat,X,G,'T',kX)
%  can also be used, with meanings as in mds_add_prod.
function FS = setprod(nparmat,X,Y,varargin)
  r = size(X,1);
  if ~isstruct(Y), n=size(Y,2); else, n=size(Y.mat,2); end
  FS = mds_set_zero(nparmat,r,n);
  FS = mds_add_prod(FS,X,Y,varargin{:});
end
