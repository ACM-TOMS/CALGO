% OMEGA_BUILD  Build band covariance matrix for a VARMA model
% 
%  [Su,Olow] = OMEGA_BUILD(S,G,W,p,n) finds cov. matrix Omega of w=(w1'...wn')':
%  
%    Su    h·r×h·r, upper left part of Omega, where h = max(p,q)
%    Olow  (n-h)r×(q+1)r, block-diagonals of lower part of Omega
%    S     r×(p+1)·r, = {S0 S1...Sp}, Si = cov(x(t),x(t-i))
%    G     r×(q+1)r, = {G0...Gq}, Gi = cov(y(t),y(t-i))
%    W     r×(q+1)r, = {W0...Wq}, Wi = cov(y(t),x(t-i))
%    n     length of series
%    p     number of autoregressive terms
%    q     number of moving average terms
%    r     dimension of xt
% 
%  Here is the structure of Omega for q<p:
% 
%    S0 S1'...      Sp-1'|              transpose
%    S1 S0             : |               of lower
%    S2 S1 S0            |              left part    r·p
%     :        ...    S1'|
%    Sp-1             S0 |
%    --------------------+-----------------------
%          Gq ...  G2 G1 | W0           transpose
%             Gq      G2 | W1  W0        of lower
%                ..   G3 | W2  W1  W0        part
%                  ..  : |  :   :    ..
%                     Gq | Wq-1        ..            r·(n-p)
%                        | Wq            ..
%                        |     Wq          ..
%                        |        ...        ..
%                        |            Wq ..... W0
%  
%  For q >= p the picture is:
%  
%    S0 S1'...   Sp-1' Gp'...  Gq-1'|                transpose
%    S1 S0         :    :         : |                 of lower
%    S2     ..     :    :         : |                left part
%     :        .. S1'  G2'        : | 
%    Sp-1      S1 S0   G1'...  Gq-p'|                            r·h
%    Gp        G2 G1   W0 ...Wq-p-1'|
%    :             :    :  ..     : |
%    :             :    :     ..  : |
%    Gq-1...     Gq-p Wq-p-1...  W0 |
%    -------------------------------+-------------------------
%         Gq .. Gq-p+1 Wq-p ...  W1 | W0             transpose
%             ..   :    :        W2 | W1  W0          of lower
%                 Gq    :         : |  :    ..            part
%                      Wq         : |  :       ..               r·(n-h)
%                           ..    : |  :          ..
%                                Wq |  :             ..
%                                   | Wq                ..
%                                   |      ..             ..
%                                   |          Wq ...       W0
%
%  The Olow matrix is r·(n-p)×r·(q+1) and the following example for q = 4, n-p =
%  5 shows its structure:
%  
%    G5 G4 G3 W2 W1 W0
%    G5 G4 W3 W2 W1 W0
%    G5 W4 W3 W2 W1 W0  r·(n-p)
%    W5 W4 W3 W2 W1 W0
%    W5 W4 W3 W2 W1 W0
%         r·(q+1)
%
function [Su, Olow] = omega_build(S,G,W,p,n)
  q = length(G) - 1;
  h = max(p,q);
  S = cell2mat(fliplr(S(1:p)));
  G = cell2mat(fliplr(G));
  W = cell2mat(fliplr(W));
  r = size(G,1);
  Su = zeros(h,h);
  Olow = zeros(r*(n-h), r*(q+1));
  K = 1:r;
  for t = 1:p
    Su(K,1:t*r) = S(:,end-t*r+1:end);
    K = K + r;
  end
  J = (q-p)*r+1 : q*r;
  for t = p+1:h
    Su(K,1:p*r) = G(:,J);
    Su(K,p*r+1:t*r) = W(:,end-(t-p)*r+1:end);
    J = J - r;
    K = K + r;
  end
  K = 1:r;
  je = min(p,q)*r;
  for t = 1:n-h
    Olow(K,1:je) = G(:,1:je);
    Olow(K,je+1:end) = W(:,je+1:end);
    if je > 0, je = je-r; end
    K = K + r;
  end
end
