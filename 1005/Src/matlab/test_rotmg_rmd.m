function test_rotmg_rmd(varargin)
  % TEST_ROTMG_RMD() checks that rotmg_rmd computes the correct adjoints. 
  % TEST_ROTMG_RMD(1) prints out a table of the rotmg flag values and those of
  %   the scalings gamma that are tested.
  %
  % NOTES: Rotmg_rmd has several if- and while-statements, and the execution
  % path depends on the values of flag, as well as the scaling done in rotmg.
  % Thus the testing of rotmg_rmd includes several tests, to ensure full code
  % coverage, and that several execution paths are tested.
  %
  % It turns out that one must be careful in computing the numerical derivatives
  % to compare with. Firstly it is necessary to choose the values where the
  % adjoints are computed so that they are not too close to points where rotmg is
  % not differentiable, secondly one must use a finite difference stencil that
  % produces accurate results: a four-point stencil transpires to be sufficient,
  % and thirdly one must choose the discretization step carefully.
  
  mini = 1e-18;
  small = 1e-9;
  big = 1e9;
  maxi = 1e17;
  testu = [
    2      2      2    3
    2      2      3    2
    small  small  2    3
    small  small  3    2
    small  mini   3    2
    mini   mini   1    3
    mini   mini   1    3
    small  2      2    2
    small  2      2e4  2
    2      small  2    2e4
    2      small  2    2
    big    big    2    3
    big    big    3    2
    maxi   big    2    3
    ];
  step = zeros(size(testu));
  step(:) = 4e-2;
  step(testu == small) = small/50;
  step(testu == mini) = mini/50;
  step(testu == big) = big/50;
  step(testu == maxi) = maxi/50;
  if nargin > 0, disp('Flag  Gamma1       Gamma2'), end
  rng(47);
  for t = 1:size(testu, 1)
    u = testu(t,:)';
    tol = 1e-5;
    for i = 1:1 % Basic testcase plus two additional random ones
      if nargin > 0, scaling_table(u), end
      ftest(@f, @f_rmd, u          ...
        ,   'tol',      tol        ...
        ,   'accum' ,   false(4,1) ...
        ,   'F2Gpar',   1          ...
        ,   'step',     step(t,:)  ...
        ,   'testcols', 1 ...
        ,   'fourpoint', [true true true true])
      u = testu(t,:)'.*(1 + (rand(4,1)-0.5)/4);
      tol = 1e-4;
    end
  end
end

function [v, flag] = f(u)
  % Translate general function evaluation with f to specific rotmg call
  d = u(1:2);
  x = u(3:4);
  [d, x1, param] = rotmg(d, x);
  flag = param(1);
  v = [d; x1];
  switch flag
    case 1, v(4:5) = param([2,5]);
    case 0, v(4:5) = param(3:4);
    case -1, v(4:7) = param(2:5);
    otherwise, assert(flag == -2) 
  end
end

function ua = f_rmd(~, v, ~, va, flag)
  % Translate general f_rmd to specific rotmg_rmd
  dt = v(1:2);
  x1t = v(3);
  dta = va(1:2);
  x1ta = va(3);
  switch flag
    case 1
      param = [flag, v(4), 0, 0, v(5)];
      parama = [flag, va(4), 0, 0, va(5)];
    case 0
      param = [flag, 0, v(4), v(5), 0];
      parama = [flag, 0, va(4), va(5), 0];
    case -1
      param = [flag, v(4:7)'];
      parama = [flag, va(4:7)'];
    case -2
      param = flag;
      parama = flag;
  end
  [da, xa] = rotmg_rmd(dt, x1t, param, dta, x1ta, parama);
  ua = [da, xa];
end

function scaling_table(u)
  % Print table of combinations of Flag and scalings that are tested
  % (used to ascertain that all code segments in rotmg_rmd are actually tested)
  d = u(1:2);
  x = u(3:4);
  [d, x1, param] = rotmg(d, x);
  dta = ones(size(d));
  x1ta = 1;
  parama = [param(1), 1, 1, 1, 1];
  [~, ~, info] = rotmg_rmd(d, x1, param, dta, x1ta, parama);
  fprintf('%-5d %-12.6g %-.6g\n', info.flag, info.gam1, info.gam2);
end
