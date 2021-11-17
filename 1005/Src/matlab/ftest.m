% FTEST  Test implementation of reverse mode derivative calculation
%
%   FTEST(@F, @G, u) compares numerical and analytic Jacobians of a function to
%   assert whether the analytic Jacobian has been implemented correctly. The
%   function F should have signature:
%
%        function v = F(u)
%
%   where u and v are (row or column) vectors. Furthermore it is assumed that
%   the setting of the computation is that a scalar f is the final result of the
%   computation and that the gradient df/dx is sought. The function G should
%   implement reverse mode differentiation of F and have signature:
%
%        function ua = G(u, v, ua, va)
%
%   where v has been calculated earlier with v = F(u), va = df/dy and df/dx is
%   added to ua (the i signifies "adjoint"). Graphically the computations can
%   be described:
%                       _________________________________
%                      /                                 \
%     u ---> F(.) ---> v --...--> f --...--> va ---> G(u,v,ua,.) ---> ua
%      \______________________________________________________/
%
%   To check the computations, we set f = yj for all j in turn and call the
%   functions with u, and also with u + d*ek and u - d*ek for all k (where d is
%   small and ek is the k-th unit vector). This will generate both a complete
%   numerical and analytic Jacobian. If the maximum difference between these
%   Jacobians, either relative or absolute, is less than a tolerence, tol, then
%   FTEST will exit quietly, otherwise an error results (the absolute difference
%   is max(|J(i,j) - Jnum(i,j)]) and in the relative difference this is divided
%   by max|J(i,j)|). Default value for tol is 1e-7.
%
%   FTEST(...'tol', tol...) sets the tolerance for the testing
%
%   FTEST(...'repeats', r...) repeats the test for an additional r-1 values of u
%   generated with rand(size(u)).
%
%   FTEST(...'testcase', {par1, par2...}...) may be used to pass extra
%   parameters to F and G, which should now have signatures v = F(u, par1,
%   par2...) and ua = G(u, v, ua, va, par1, par2...).
%
%   TEST(...'accum', accum...) may be used to omit the test of accumulation for
%   some elements of ua (the ones corresponding to false elements in the logical
%   vector accum which must have same size as u). This is apropriate for testing
%   F operations that overwrite some input arguments.
%
%   TEST(...'F2Gpar', k...) causes calls to F contain k extra parameters,
%   which are passed to subsequent G-calls as the last parameters. Default k=0.
%
%   TEST(...'step', step...) sets the finite difference step for u(i) to step(i)
%   (default is 1e-5)
%
%   TEST(...'testcols', [j1 j2...]...) tests only adjoints for outputs j1, j2...
%
%   TEST(...'fourpoint', [fp1, fp2...]...) uses a four-point finite difference
%   stencil for the ua corresponding to true fpi. The default is to use a
%   two-point stencil.

function ftest(F, G, u, varargin)
  ischr = @(s) ischar(s) || isstring(s);
  p = inputParser;
  p.addParameter('tol', 1e-5)
  p.addParameter('repeats', 1)
  p.addParameter('accum', true(size(u)))
  p.addParameter('F2Gpar', 0)
  p.addParameter('testcase', {})
  p.addParameter('step', 1e-5 * ones(size(u)))
  p.addParameter('testcols', "all")
  p.addParameter('fourpoint', false(size(u)))
  p.parse(varargin{:})
  testcase = p.Results.testcase;
  assert(iscolumn(u))
  rep = 1;
  m = length(u);
  nextra = p.Results.F2Gpar;
  accum = p.Results.accum;
  fourpoint = p.Results.fourpoint;
  if isscalar(fourpoint), fourpoint = repmat(fourpoint, size(u)); end
  
  % RUN TESTS rep TIMES:
  while rep <= p.Results.repeats
    
    % CHECK JACOBIAN:
    [v, extra{1:nextra}] = F(u, testcase{:});
    n = length(v); % the Jacobian will be m by n with J(i,j) = dxi/dyj
    jlist = p.Results.testcols;
    if ischr(jlist) && jlist == "all", jlist = 1:n; end
    for i = 1:m
      d = p.Results.step(i);
      Fm = F(u - d*e(i), testcase{:});
      Fp = F(u + d*e(i), testcase{:});
      if fourpoint(i)
        Fmm = F(u - 2*d*e(i), testcase{:});
        Fpp = F(u + 2*d*e(i), testcase{:});
        Jnum(i,:) = (Fmm(:) - Fpp(:) + 8*(Fp(:) - Fm(:)))'/(12*d);
      else
        Jnum(i,:) = (Fp(:) - Fm(:))'/(2*d);
      end
    end
    k = 0;
    for j = jlist
      va = zeros(n, 1);
      va(j) = 1;
      ua = zeros(size(u));
      ua = G(u, v, ua, va, testcase{:}, extra{1:nextra});
      k = k+1;
      Jrmd(:,k) = ua(:);
      if any(accum)
        % CHECK THAT G ACCUMULATES IN ua:
        ui1 = zeros(size(u));
        ui1(p.Results.accum) = 1;
        ui1K = ui1;
        ui1 = G(u, v, ui1, va, testcase{:}, extra{1:nextra});
        OK = norm(ui1 - ui1K - ua) < 1e-10;
        assert(OK, 'G does not accumulate correctly');
      end
    end
    Jnum = Jnum(:,jlist);
    assert(isequal(size(Jnum),size(Jrmd)));
    Jdiff = Jrmd - Jnum;
    reldiff = abs(Jdiff)./max(1e-12,abs(Jrmd));
    if max(reldiff(:)) > p.Results.tol
      disp('Analytic Jacobian'), disp(Jrmd)
      disp('Numerical Jacobian'), disp(Jnum)
      disp('Relative difference'), disp(reldiff)
      disp('Max relative difference'), disp(max(reldiff(:)))
      error('Jacobians do not match');
    end
    %disp('Relative difference'), disp(reldiff)
    %disp('Max relative difference'), disp(max(reldiff(:)))
  
    % CHECK THAT G CAN HANDLE MULTIPLE FINAL OUTPUTS
    %     nf = 1;
    %     a = randPM(nf,n);
    %     %f = a*v;
    %     va = a';
    %     ua = zeros(m, nf);
    %     ua = G(u, v, ua, va, varargin{:});
    %     ui1 = J*a';  % m x nf
    %     if ~almostequal(ua, ui1)
    %       error('Check of multiple final output fails');
    %     end
    
    % PREPARE FOR NEXT REPEAT
    rep = rep + 1;
    u = rand(size(u));
  end
  
  function e = e(i)
    e = zeros(m, 1);
    e(i) = 1;
  end
  
end
