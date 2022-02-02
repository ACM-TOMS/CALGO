% TEST_VYW  Tests solution of vector-Yule-Walker system
%
%   TEST_VYW QUIET prints out "Passed vyw_tests OK" if all tests are passed.
%   TEST_VYW (without an argument) prints more information.
%   TEST_VYW(i:j) runs testcases number i through j.
%   TEST_VYW MRHS tests multiple rhs feature. 
%   TEST_VYW MRHS QUIET, TEST_VYW('MRHS',i) etc. may also be used.
%
%   See also: testcase
%
function test_vyw(varargin)
  [varargin,quiet,mrhs] = getflags(varargin,'quiet','mrhs');
  fprintf('TESTING VYW_FACTORIZE AND VYW_SOLVE...'); 
  fprintf_if(~quiet,'\n');
  A={}; B={}; Sig={}; name={};
  if ~isempty(varargin) && isnumeric(varargin{1})
    cases = varargin{1};
    for i=1:length(cases)
      [A{end+1},B{end+1},Sig{end+1},name{end+1}] = testcase(cases(i));
    end
    nNamed = 0;
  else
    [A,B,Sig,name] = testcase('all'); % OBTAIN NAMED TESTCASES
    nNamed = length(A);
    for p=0:3 % ADD UNNAMED CASES FOR SEVERAL COMBINATIONS OF p,q,r
      for q=0:3
        for r=1:3
          [A{end+1},B{end+1},Sig{end+1}] = testcase(p,q,r);
          name{end+1}='Unnamed';
        end
      end
    end
  end
  for i = 1:length(A)
    % PRINT INFO
    r = size(Sig{i},1);
    p = size(A{i},2)/r; 
    q = size(B{i},2)/r; 
    fprintf_if(~quiet,'\n  %s %s, p=%d, q=%d, r=%d\n  ','Testcase:',name{i}...
      ,               p,q,r);
    if is_stationary(A{i}, Sig{i})
      fprintf_if(~quiet, 'stationary model\n')
    else
      fprintf_if(~quiet, 'non-stationary model\n')
    end
    [C,G,W] = find_CGW(A{i},B{i},Sig{i});
    PLU = vyw_factorize(A{i});
    S = vyw_solve(A{i},PLU,G);
    check_solution(quiet,A{i},G,S);
    if i<=nNamed && mrhs  
      % CHECK MULTIPLE RHS FEATURE FOR NAMED TESTCASES:
      Bc = cell(1,q); Gc = cell(1,q+1);
      for c=1:3
        for j=1:q, Bc{j} = B{i}{j}-0.1*c; end
        Sigc = Sig{i} + 0.1*c;
        [C,Gc,W] = find_CGW(A{i},Bc,Sigc);
        for j=1:q+1, G{j}(:,:,c) = Gc{j}; end
      end
      PLU = vyw_factorize(A{i});  % SOLVE
      S = vyw_solve(A,PLU,G);  % SOLVE
      fprintf_if(~quiet,'MRHS: ')
      check_solution(quiet,A{i},G,S);
    end
  end
  disp('OK')
end

function check_solution(quiet,A,G,S)
  % Check residuals of original equations
  q = length(G) - 1;
  r = size(G{1},1);
  p = size(A,2)/r;
  N = size(G{1},3);
  maxres = 0;
  Sc = {};
  for c=1:N
    for i=1:p+1, Sc{i} = S{i}(:,:,c); end
    for i=1:q+1, Gc{i} = G{i}(:,:,c); end
    maxres = max(maxres, find_max_residual(A,Sc,Gc));
  end
  fprintf_if(~quiet,'Max residual = %.1e',maxres)
  ascertain(maxres < 1e-12);
  fprintf_if(~quiet,'...OK\n  ');
end

function maxres = find_max_residual(A,S,G)
  % Max residual of modif. VYW-system. A=[A1...Ap], S={S0...Sp}, G={G0...Gq}
  q = length(G) - 1;
  r = size(G{1},1);
  A = makecell(A);
  p = length(A);
  S0 = S{1}; S = S(2:end);
  G0 = G{1}; G = G(2:end);
  maxres = 0;
  for i = 0:p
    if     i == 0, sum = S0 - G0  ;
    elseif i <= q, sum = S{i} - G{i};
    else         , sum = S{i}       ;
    end
    for j = 1:p
      if j <  i, sum = sum - A{j}*S{i-j} ; end
      if j == i, sum = sum - A{j}*S0      ; end
      if j >  i, sum = sum - A{j}*S{j-i}'; end
    end
    maxres = max(maxres,norm(sum));
  end
end
