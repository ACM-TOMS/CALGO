% ATB_C  Calculate A'·B + C and optionally its derivatives
%
%  D = ATB_C(A,B,C) calculates D = A'·B + C.
%
%  [D,Dd] = ATB_C(A,B,C,Ad,Bd,Cd) calculates also the derivative w.r.t.
%  parameter theta(i) in Dd(:,:,i), i = 1,...,nPar. Ad(:,:,i), Bd(:,:,i) and
%  Cd(:,:,i) should be the derivatives of A, B and C w.r.t. theta(i).
%
%  Dd = ATB_C(A,B,Ad,Bd,Cd) calculates only the derivatives.
%
%  When C is zero use ATB_C(A,B), and ATB_C(A,B,Ad,Bd)
%
%  METHOD: Dd(:,:,i) = Ad(:,:,i)'*B + A'*Bd(:,:,i) + Cd.

function [varargout] = atb_c(A, B, varargin)  
  DIFF = nargin > 3;
  if DIFF
    switch nargin
      case 4, [Ad,Bd] = deal(varargin{:}); C = []; Cd = [];
      case 5, [Ad,Bd,Cd] = deal(varargin{:}); C = [];
      case 6, [C,Ad,Bd,Cd] = deal(varargin{:});
      otherwise error('illegal call');
    end
    [k,m] = size(A);
    n = size(B,2);
    nPar = size(Ad,3);
    Dd = reshape(A'*reshape(Bd,k,n*nPar), m, n, nPar) + ...
      permute(reshape(B'*reshape(Ad,k,m*nPar), n, m, nPar),[2,1,3]);
    if ~isempty(Cd), Dd = Dd + Cd; end
    varargout{nargout} = Dd;
  end
  if ~DIFF || nargout==2
    if nargin==3 || nargin==6, C = varargin{1}; else C = []; end
    D = A'*B;
    if ~isempty(C), D = D + C; end
    varargout{1} = D;
  end
end
