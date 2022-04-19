function demo_rot
  %  Demonstrate computation of adjoints of A and b for the norm of the 
  %  solution to a linear least squares problem:
  %
  %                          min ||Ax - b||
  %                           x
  %
  %  If A is QR-factorized: A = Q*R = [Q1 Q2]*[R1;0] then f is given by
  %
  %                        f = ||x||, with:
  %                        c = Q1'*b
  %                        x = R1^(-1)*c
  %
  %  If A is n+1 by n with zeros below the subdiagonal then R1 and Q1'b may
  %  be computed by n Givens rotations in the (1,2),...(n,n+1) planes that
  %  zero the (2,1)...(n+1,n)-elements in the matrix H = [A b]: After the
  %  application of the rotations H will contain [R Q'*b]. A subsequent
  %  trsv-operation provides x.
  %
  %          H = | x  x  x  x |  Q'H = | x  x  x  x | = [R c]
  %              | x  x  x  x |        |    x  x  x |
  %              |    x  x  x |        |       x  x |
  %              |       x  x |        |          x |
  %
  %  x and f are examples of functions of Givens rotations that are
  %  continuously differentiable, despite the fact that the rotations
  %  themselves are not continuous.
  
  A = [ 
    1  2
   -3  1
    0  4
    ];
  b = [5 8 8]';
  [f, Aa, ba] = example2(A, b);
  fprintf('\nf = %.6f\n', f); 
  print_adjoints('Analytic adjoints', Aa, ba)
  test_example2_rmd(A, b); % Compare with numerical derivatives
end

function [f, Aa, ba] = example2(A, b)
  % Juggle arguments (for use with fnumd)
  a = [A b];
  
  % Zero h21. Row 2 rotated is kept in z for use in the rmd operation. In the
  % general n+1 by n case then n-1 such storage vectors are needed.
  [G1, a(1,1)] = rotg(a(1,1), a(2,1));
  [a(1,2:3), a(2,2:3)] = rot(G1, a(1,2:3), a(2,2:3));
  z = a(2,2:3); % Store for rmd operation

  % Zero h32
  [G2, a(2,2)] = rotg(a(2,2), a(3,2));
  [a(2,3), a(3,3)] = rot(G2, a(2,3), a(3,3));
  
  % Compute function value
  R1 = triu(a(1:2,1:2));
  c = a(1:2,3);
  x = R1\c;
  f = norm(x);
  
  % RMD COMPUTATION:
  if nargout > 1
    % Rmd for norm(x)
    fa = 1;
    xa = fa*x/f;
    
    % Rmd for trsv (x = R1\c)
    ca = R1'\xa;
    R1a = -triu(ca*x');
    
    % Store ca and R1a in ha
    ha = zeros(size(a));
    ha(1:2,3) = ca;
    ha(1:2, 1:2) = R1a;
    
    % Rmd for (2,3) plane rotation
    [ha(2,3), ha(3,3), c2a, s2a] = rot_rmd(a(2,3), a(3,3), G2, ha(2,3),ha(3,3));
    [ha(2,2), ha(3,2)] = rotg_rmd(G2, a(2,2), c2a, s2a, ha(2,2));
    
    % Rmd for (1,2) plane rotation
    I = 2:3;
    [ha(1,I), ha(2,I), c1a,s1a] = rot_rmd(a(1,I), z, G1, ha(1,I), ha(2,I));
    [ha(1,1), ha(2,1)] = rotg_rmd(G1, a(1,1), c1a, s1a, ha(1,1));
    
    % Extract Aa and ba
    [Aa, ba] = deal(ha(:,1:2), ha(:,3));
  end
end

% FUNCTION FOR CHECKING THAT THE ADJOINTS IN EXAMPLE 2 ARE CORRECT
function test_example2_rmd(A, b)
  repeats = 100;
  for rep = 1:repeats
    x = [A(1:2,1); A(1:3,2); b(:)];
    n = length(x);
    
    % Compute numerical adjoints with central finite differences:
    a = 1e-6;
    for i = 1:n
      fa(i,1) = (f(x + a*e(i)) - f(x - a*e(i)))/(2*a);
    end
    [normx, Aa, ba] = example2(A, b);
    assert(almostequal(normx, norm(A\b))); % Compare with Matlab's LS solution

    fi_analytic = [Aa(1:2,1); Aa(1:3,2); ba(:)];
    
    % Print numerical adjoints on first round:
    if rep == 1
      numAa(1:2,1) = fa(1:2);
      numAa(1:3,2) = fa(3:5);
      print_adjoints('Numerical adjoints', numAa, fa(6:8))
    end
    if max(abs(fa - fi_analytic)) > 0.1, keyboard, end
    assert(almostequal(fa, fi_analytic, 1e-4))
    A = 10*rand(size(A)); b = 10*rand(size(b)); % for next loop
    A(3,1) = 0;
    if cond(A) > 50, continue, end
  end
  function f = f(x)
    A = [[x(1:2);0], x(3:5)];
    b = x(6:8);
    f = example2(A, b);
  end
  function e = e(i), e = zeros(n,1); e(i) = 1; end
end

function print_adjoints(s, Aa, ba, dig)
  if nargin < 4, dig = 4; end
  fprintf('\n%s:\nAa = \n', s);
  fm = sprintf('%%%d.%df', dig+3, dig);
  fm1 = [fm ' ' fm '\n'];
  fm2 = ['\nba =\n' fm '\n' fm '\n' fm '\n'];
  for i=1:3, fprintf(fm1, Aa(i,:)); end  
  fprintf(fm2, ba);
end

% MATLAB IMPLEMENTATION OF THE BLAS ROTG AND ROT  
function [x, y] = rot(G, x, y)
  [c, s] = deal(G(1,1), G(1,2));
  [x,y] = deal(c*x + s*y, -s*x + c*y);
end

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
