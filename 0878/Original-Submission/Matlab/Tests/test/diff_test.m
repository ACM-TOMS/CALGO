%DIFF_TEST  Check analytic derivative against numerical derivative
%
%  DMAX = DIFF_TEST(@FUN, X, PAR1, PAR2,...) checks the derivative of FUN at X.
%  FUN should be a function with header [F, G] = FUN(X, PAR1, PAR2,...) where
%  X is an m-vector, F is an n-dimensional function value and G is an n×m
%  matrix with the Jacobian of FUN at X. DMAX returns the maximum relative
%  difference between G(i,j) and a numerical derivative calculated with the
%  central difference formula.
%
%  [DMAX,G,GNUM] = DIFF_TEST(...) returns also the G from FUN and the numerical
%  Jacobian.
%
%  For efficiency it is possible to test a FUN that returns several functions
%  and gradients in cell arrays F = {f1...fn} and G = {g1...gn} where fi and gi
%  are the function value and Jacobian of the i-th function. With this syntax
%  DMAX is a vector of the discrepancies.

function [dmax,g,gnum] = diff_test(fun, x, varargin)
  [f,g] = fun(x, varargin{:});
  if ~iscell(f), f={f}; g={g}; end
  Nf = length(f);
  n = length(x);
  eps = norm(x)*1e-6;
  xp = x;
  xm = x;
  for i=1:n
    xp(i) = x(i)+eps;
    xm(i) = x(i)-eps;
    fp = fun(xp, varargin{:});
    fm = fun(xm, varargin{:});
    if ~iscell(fp), fp = {fp}; fm = {fm}; end
    for j=1:Nf
      gnum{j}(:,i) = (fp{j} - fm{j})/(2*eps);
    end
    xp(i) = x(i);
    xm(i) = x(i);
  end
  for j=1:Nf
    ascertain(isequal(size(g{j}), size(gnum{j})));
    ascertain(all(isfinite(g{j}(:))));
    ascertain(all(isfinite(gnum{j}(:))));
    gmax = max(max(abs(g{j})));
    if isempty(gnum{j}) && isempty(g{j})
      dmax(j) = 0;
    elseif gmax<1e-13
      dmax(j) = max(max(abs(gnum{j}-g{j})));
    else
      dmax(j) = max(max(abs(gnum{j}-g{j})))/gmax;
    end
  end 
end
