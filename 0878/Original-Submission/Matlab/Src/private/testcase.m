% TESTCASE  Create testcases for VARMA likelihood calculation
%
%   [A,B,Sig] = TESTCASE(NAME) returns a named test case in matrices A =
%   [A1...Ap], B = [B1...Bq] and Sig suitable for testing various components of
%   the calculation of exact VARMA likelihood. TESTCASE SUMMARY prints a list of
%   the named testcases (with possible values for NAME)
%
%   [A,B,Sig,p,q,r] = TESTCASE(NAME) returns also the dimensions of the case.
%
%   [A,B,Sig,name] = TESTCASE(i) returns the i-th named testcase. To get also
%   dimensions use [A,B,Sig,p,q,r,name] = TESTCASE(i).
%
%   n = TESTCASE('NUMBER') returns the number of different named testcases.
%
%   [A,B,Sig,name] = TESTCASE('ALL') returns all the named testcases in three
%   cell arrays; the i-th one is returned in A{i}, B{i} and Sig{i}. In addition
%   name{i} returns the name of the i-th case.
%
%   [A,B,Sig] = TESTCASE(p,q,r) returns an unnamed testcase with dimensions
%   p,q,r.

function [A,B,Sig,varargout] = testcase(varargin)
  name = '';
  if nargin == 1
    type = varargin{1};
  else
    p = varargin{1};
    q = varargin{2}; 
    r = varargin{3};
    type = 'unnamed';
  end
  testcases = {
    'tinyAR'
    'tinyMA'
    'tinyARMA'
    'smallAR'
    'smallMA'
    'smallARMA1'
    'smallARMA2'
    'mediumAR'
    'mediumARMA1'
    'mediumARMA2'
    'mediumMA'
    'largeAR'
    %'pivotfailure'
    };
  N = length(testcases);
  if ~ischar(type), name = testcases{type}; 
  elseif ~isequal(type,'all'), name=type; end
  switch type
    case 'number'
      A = N; % count of named cases
      return
    case 'summary'
      fprintf('No. Name          p  q  r  Stationary?\n');
      for i = 1:N
        fmt = '%-2i  %-12s %2d %2d %2d   %6s\n';
        [A,B,Sig,name] = testcase(i);
        [r,p,q] = get_dimensions(A,B,Sig);
        if is_stationary(A, Sig) stat = 'yes'; else stat = 'no'; end
        fprintf(fmt, i, name, r, p, q, stat);
      end
      clear A B Sig name
      return
    case 'all'
      for i=1:N
        name{i} = testcases{i};
        [A{i}, B{i}, Sig{i}, p{i}, q{i}, r{i}] = testcase(name{i});
      end
      if any([4,7] == nargout), varargout{nargout-3} = name; end
      if nargout > 4, [varargout{1:3}] = deal(p,q,r); end
      return
    case 'unnamed'
      A = rand(r,r*p)/(p*r)*2;
      B = rand(r,r*q)/(q*r)*2;
      Sig=hilb(r)+0.2*eye(r);
      while ~is_stationary(A,Sig), A=A/2; end
    case {1,'tinyAR'}
      A = 0.5;
      B = [];
      Sig = 0.8;
    case {2,'tinyMA'}
      A = [];
      B = 0.5;
      Sig = 0.8;
    case {3,'tinyARMA'}
      A = 0.4;
      B = 0.4;
      Sig = 0.8;
    case {4,'smallAR'}
      %A = [[0.3 0.1; 0.4 0.2],[0.1 0.1;0.1 0.2]];
      A = [0.1 0.1; 0.1 0.1];
      B = [];
      Sig = [2 1;1 3];
    case {5,'smallMA'}
      A = [];
      B = [0.3 0.1 0.1 0.1; 0.3 0.1 0.1 0.1];
      Sig = [2 1;1 2];
    case {6,'smallARMA1'}
      A = [0.3 0.1; 0.4 0.2];
      B = [0.2 0.3; 0.2 0.3];
      Sig = [2 1;1 2];
    case {7,'smallARMA2'}
      B = [[0.4 0.1; 0.2 0.1], [0.2 0.1; 0.3 0.1]];
      A = [0.2 0.3; 0.2 0.3];
      Sig = [2 1;1 2];
    case {8,'mediumAR'}
      A = [
        0.35 0.15 0.15  0.11 0.14 0.07;
        0.25 0.15 0.05  0.12 0.15 0.08;
        0.15 0.05 0.01  0.13 0.16 0.09];
      B = [];
      Sig=[2.0  0.5 0.0;
           0.5  2.0 0.5;
           0.0  0.5 1.0];
    case {9,'mediumARMA1'}
      A = [
        0.15 0.10 0.05  0.11 0.14 0.17  0.01 0.04 0.06;
        0.16 0.11 0.06  0.12 0.15 0.18  0.02 0.05 0.08;
        0.17 0.12 0.07  0.13 0.16 0.19  0.03 0.06 0.09];
      B = fliplr(A);
      Sig=[2.0  0.5 0.0;
           0.5  2.0 0.5;
           0.0  0.5 1.0];
    case {10,'mediumARMA2'}
      B = [
        0.15 0.10 0.05  0.11 0.14 0.17  0.01 0.04 0.06;
        0.16 0.11 0.06  0.12 0.15 0.18  0.02 0.05 0.08;
        0.17 0.12 0.07  0.13 0.16 0.19  0.03 0.06 0.09];
      A = fliplr(B);
      Sig=[2.0  0.5 0.0;
           0.5  2.0 0.5;
           0.0  0.5 1.0];
    case {11,'mediumMA'}
      A = [];
      B = [
        0.35 0.25 0.15  0.11 0.14 0.17;
        0.25 0.15 0.05  0.12 0.15 0.18;
        0.15 0.05 0.01  0.13 0.16 0.19];
      Sig=[2.0  0.5 0.0;
           0.5  2.0 0.5;
           0.0  0.5 1.0];
    case {12,'largeAR'}
      r=7; p=5; rand('state',6);
      A = 1.8*rand(r,r*p)/r/p;
      Sig=hilb(r)+0.2*eye(r);
      B = [];
    case {'pivotfailure'} % create almost singular vyw equations
      A = [
        -0.6250  0.2400    0.2500    0.1250    0.6250    0.3750
        0             0    0.1250    0.5000    0.2500    0.1250
        0.2500        0    0.1430         0    0.2500    0.2500];
      A(1,1) = -0.62999813363195223; %% as singular as can be made
      B=[];
      Sig = eye(3);
    otherwise
      error('no such testcase')
  end
  if any([4,7] == nargout), varargout{nargout-3} = name; end
  if nargout > 4
    r = size(Sig,1);
    p = size(A,2)/r;
    q = size(B,2)/r;
    [varargout{1:3}] = deal(p, q, r);end
  ascertain(isequal(Sig,Sig'));
  ascertain(min(eig(Sig))>0);
end
