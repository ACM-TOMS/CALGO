function test_rot
  % Check that formulae for the adjoints of Givens rotation generation and
  % application (ROT and ROTG) are correct,
  test_rmd()
end

function test_rmd
  for n = 1:4
    c = rand();
    s = sqrt(1-c^2);
    x = rand(2*n, 1);
    u = [c; s; x];
    accum = [true, true, false(1, 2*n)];
    ftest(@frot, @frot_rmd, u, 'accum', accum)
    u = rand(2,1);
    ftest(@frotg, @frotg_rmd, u, 'accum', [false, false])
  end
end

function v = frotg(u)
  a = u(1);
  b = u(2);
  [G, x1] = rotg(a, b);
  v(1) = G(1,1);
  v(2) = G(1,2);
  v(3) = x1;
end

function ua = frotg_rmd(~, v, ~, va)
  [c, s] = deal(v(1), v(2));
  [ca, sa] = deal(va(1), va(2));
  x1 = v(3);
  x1a = va(3);
  G = [c, s; -s, c];
  [aa, ba] = rotg_rmd(G, x1, ca, sa, x1a);
  ua = [aa;ba];
end

function v = frot(u)
  c = u(1);
  s = u(2);
  n = (length(u) - 2)/2;
  x = u(3:n+2)';
  y = u(n+3:end)';
  G = [c, s; -s, c];
  [x, y] = rot(G, x, y);
  v = [x(:); y(:)];
end

function ua = frot_rmd(u, v, ua, va)
  n = (length(u) - 2)/2;
  [c, s] = deal(u(1), u(2));
  G = [c, s; -s, c];
  [ca, sa] = deal(ua(1), ua(2));
  x = v(1:n)';
  y = v(n+1:end)';
  xa = va(1:n)';
  ya = va(n+1:end)';
  [xa, ya, ca, sa] = rot_rmd(x, y, G, xa, ya, ca, sa);
  ua = [ca; sa; xa(:); ya(:)];
end
  
% MATLAB IMPLEMENTATION OF THE BLAS ROTG AND ROT  
function [x, y] = rot(G, x, y)
  [c, s] = deal(G(1,1), G(1,2));
  [x,y] = deal(c*x + s*y, -s*x + c*y);
end
%
function [G, x1] = rotg(a,b)
  % Translation of the reference BLAS srotg
  % Return G s.t. G*[a;b] = [x1;0]
  r = hypot(a, b);
  if r == 0
    c = 1;
    s = 0;
  else
    if abs(a) > abs(b) && a < 0 || abs(a) <= abs(b) && b < 0
      r = -r;
    end
    c = a/r;
    s = b/r;
  end
  G = [c s; -s c];
  x1 = r;
end

% MATLAB IMPLEMENTATION OF REVERSE MODE DERIVATIVES FOR ROTG AND ROT
function [xa,ya] = rotg_rmd(G, x, ca, sa, xa)
  [c,s] = deal(G(1,1), G(1,2));
  if x ~= 0
    c1 = ca/x;
    s1 = sa/x;
    d1 = xa - s*s1 - c*c1;
    xa = c1 + c*d1;
    ya = s1 + s*d1;
  end
end
%
function [xa, ya, ca, sa] = rot_rmd(x, y, G, xa, ya, ca, sa)
  if nargin < 7, ca = 0; sa = 0; end
  % x and y are values on exit from rot; row vectors
  [c, s] = deal(G(1,1), G(1,2));
  c1 = xa*x' + ya*y';
  s1 = xa*y' - ya*x';
  ca = ca + c1*c - s1*s;
  sa = sa + s1*c + c1*s;
  [xa, ya] = deal(c*xa - s*ya, s*xa + c*ya);
end
