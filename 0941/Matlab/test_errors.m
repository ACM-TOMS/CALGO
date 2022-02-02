function test_errors(disp_message)
%TEST_ERRORS Test error handling
%
% TEST_ERRORS() calls most functions of the toolbox, and is
% intended to trigger all error(...) commands in the toolbox.
% This function is called by TEST_FUNCTIONS.
%
% TEST_ERRORS(TRUE) additionally displays all generated error messages.
%
% See also TEST_FUNCTIONS, CHECK_FUNCTIONS.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin == 0 || ~islogical(disp_message) || ~isscalar(disp_message))
  disp_message = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\ndematricize:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  dematricize(1, 2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  dematricize('a', 2, 3);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  dematricize(1, -1, 4);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  dematricize(1, 2, [1 1]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  dematricize(1, [2 2 2], [1 2], [1 3]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  dematricize(1, [2 2 2], [1 2], [1 1]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\ndiag3d:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  diag3d([1 2; 3 4]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\ngen_invlaplace:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  gen_invlaplace(1, 2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_invlaplace({[1 2]}, 2, 4);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_invlaplace(2, 2, 4);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\ngen_laplace:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  gen_laplace(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_laplace(1, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_laplace(3, {1});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_laplace(1, {[1,2]});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_laplace(3, [1,2]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_laplace(3, 'foo');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_laplace(3, 1, {1});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_laplace(1, 1, {[1,2]});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_laplace(3, 1, [1,2]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_laplace(3, 1, 'foo');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\ngen_sin_cos:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  gen_sin_cos();
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  gen_sin_cos(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nhtenones:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htenones();
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htenones(-2.2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htenones([]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nhtenrandn:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htenrandn();
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htenrandn([1 2; 3 4]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htenrandn([]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htenrandn([1 2], 3);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htenrandn([1 2], 'o', [3 4 5 6]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nkhatrirao_aux:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  khatrirao_aux(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  khatrirao_aux('a', 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  khatrirao_aux(randn(2, 2), randn(2, 3));
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nkhatrirao_t:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  khatrirao_t(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  khatrirao_t('a', 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  khatrirao_t(randn(2, 2), randn(3, 2));
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nlaplace_core:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  laplace_core(1.1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nmatricize:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  matricize(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  matricize('a', 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  matricize(1, -1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  matricize(1, 1, -2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  matricize(randn(2, 3, 4), [1 2], [3 1]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nreciproc_sum:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  reciproc_sum(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  reciproc_sum(-1, 2, 3, 4);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nspy3:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  spy3('a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\ntruncate:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  truncate(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  truncate(1, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  truncate('a', struct('max_rank', 50));
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nttm:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  ttm(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm('a', 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(1, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(1, {'a'});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(1, 2, 0.5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(1, 2, -5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(1, 2, [1 2]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(1, [2 3], 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(1, [2; 3], 1, 't');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(1, [2; 3], 1, 'h');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\nttt:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  ttt(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt('a', 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt(1, 2, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt(1, [2 3], 1, 2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/htensor:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor({});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor({'a'});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 3 4 5], [1 2 3], [4 5 6 7]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 3 4 5], [2 3; 4 5; 6 7; 0 0; 0 0; 0 0; 0 0], [4 5]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 3 -4 5]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor({randn(4, 3), randn(5, 2)});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor({randn(4, 3), randn(5, 3)}, [2 3; 4 5; 6 7; 0 0; 0 0; 0 0; 0 0], [4 5]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor({randn(4, 3), randn(5, 3)}, [2 3; 0 0; 0 0], [4 3 2 1]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor('a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/check_htensor:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


std_children = [2 3; 4 5; 6 7; 0 0; 0 0; 0 0; 0 0];
std_dim2ind = [4 5 6 7];
std_U = {[], [], [], randn(4, 3), randn(5, 3), randn(6, 3), randn(7, 3)};
std_B = {randn(4, 4), randn(3, 3, 4), randn(3, 3, 4), [], [], [], []};

try
  htensor([1 2 3], std_dim2ind, std_U, std_B);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 3; 0 0; 6 7; 0 0; 0 0; 0 0; 0 0], std_dim2ind, std_U, std_B);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 3; 4 5; 6 7; 0 0; 0 0; 0 0; 0 2], std_dim2ind, std_U, std_B);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 6; 4 5; 0 0; 0 0; 0 0; 3 7; 0 0], [4 5 3 7], std_U([1 2 6 4 5 3 7]), std_B([1 2 6 4 5 3 7]));
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, std_dim2ind', std_U, std_B);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 7 8 9], std_U, std_B);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 3], std_U, std_B);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 7], 'a', std_B);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 7], std_U', std_B);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 7], {[], [], [], randn(4, 3), randn(5, 3), randn(6, 3), randn(7, 3, 2)}, std_B);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 7], std_U, {randn(4, 4), randn(3, 3, 4), randn(3, 3, 4, 5), [], [], [], []});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 7], std_U, {randn(4, 4), randn(3, 3, 4), randn(3, 3, 5), [], [], [], []});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 7], std_U, {randn(0, 0), randn(3, 3, 0), randn(3, 3, 0), [], [], [], []});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 7], std_U, {randn(4, 4, 2), randn(3, 3, 4), randn(3, 3, 4), [], [], [], []});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor(std_children, [4 5 6 7], std_U, std_B, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/define_tree:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.define_tree();
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.define_tree([1 1]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.define_tree([1 2], 3);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.define_tree([]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.define_tree([1 3]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/subtree:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.subtree(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.subtree([1 2 3], 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.subtree([1 2; 1 2], -1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.subtree([2 1; 1 2], 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/transpose:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = htensor([3 4 5 6]);

try
  x.';
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/ctranspose:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  x';
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/disp:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  disp(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  disp(x, 'a', [1 2]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/disp_all:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  disp_all('a', x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/disp:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  isequal(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  isequal(x, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/rank:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  rank(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/size:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  size(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  size(x, 0.5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/spy:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  spy(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  spy(x, 'asdf');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/plot_sv:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  plot_sv(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  plot_sv(x, 'asdf');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/squeeze:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  squeeze(x, 0.5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/mtimes:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  'a' * x;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x * 'a';
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x * x;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/mtimes:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  'a' / x;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x / [1 2];
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/equal_dimtree:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  equal_dimtree(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  equal_dimtree(x, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/end:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  x(3, end, 4);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/permute:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  permute(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  permute(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  permute(x, 0.5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  permute(x, [1 2 3 5]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  permute(x, [1 2 1 2]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/ipermute:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  ipermute(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ipermute(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ipermute(x, 0.5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ipermute(x, [1 2 3 5]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ipermute(x, [1 2 1 2]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/subsasgn:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  x.A = 2;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x(1, 2, 3, 4) = 2;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x{0} = 2;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x.U{1} = 2;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x.B{4} = 2;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x.B{1} = randn(2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/subsref:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  x.A;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x(0.5, 2, 3, 4);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x(10, 1, 1, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x(1, 1, 1, 1, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x{1, 2};
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x{10};
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/change_dimtree:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  change_dimtree(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  change_dimtree(x, [2 3; 4 5; 0 0; 0 0; 0 0]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  change_dimtree(x, std_children, [2 3 6 7]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

x1 = htensor([4 5 6 7 8 9]);
x2 = htensor([4 5 6 7 8 9], 'TT');
try
  change_dimtree(x1, x2.children, x2.dim2ind);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/change_root:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  change_root(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  change_root(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  change_root(x, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  change_root(x, 3, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/plus:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  plus(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x + 1;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 3 4 5]) + htensor([2 3 4 5], 'TT');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 3 4]) + htensor([2 3 5]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/minus:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  minus(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x - 1;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 3 4 5]) - htensor([2 3 4 5], 'TT');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([2 3 4]) - htensor([2 3 5]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/truncate_ltr:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std_xfull = randn(3, 3, 3);
std_opts.max_rank = 5;

opts1 = std_opts;
opts1.sv = 'foo';

try
  htensor.truncate_ltr(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_ltr('a', std_opts);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_ltr(std_xfull, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_ltr(std_xfull, opts1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/truncate_rtl:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.truncate_rtl(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_rtl('a', std_opts);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_rtl(std_xfull, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_rtl(std_xfull, opts1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/truncate_std:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  truncate_std(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  truncate_std(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  truncate_std(x, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/truncate_nonorthog:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  truncate_nonorthog(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  truncate_nonorthog(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  truncate_nonorthog(x, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/truncate_cp:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.truncate_cp(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_cp(1, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_cp({'a', 'b'}, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_cp({randn(5, 2), randn(8, 2)}, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/truncate_sum:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.truncate_sum(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_sum('a', 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_sum({'a', 'b'; 'c', 'd'}, 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_sum({x1, x2}, 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.truncate_sum({x, x}, 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/trunc_rank:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.trunc_rank(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.trunc_rank('a', 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.trunc_rank(randn(4, 1), 'b');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/ttm:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  ttm(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, {'a', 'b'});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, {2, 3}, x1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, {2, 3, 4});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, {2, 3, 4}, -2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, {2, 3, 4, 5}, -8);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, {2, 3}, [2, 3, 4]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, randn(4, 4), 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, randn(4, 4), 1, 't');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttm(x, randn(4, 4), 1, 'h');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/ttt:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  ttt(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt(x, 2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt(x, x, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt(x, x, 1, 2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt(htensor([4 4 4]), htensor([4 4 4 4]), [1 2], [1 4]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt(htensor([4 4 4]), htensor([4 4 4 4 4]), [1 2 3], [1 3 5]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt(htensor([4 4 4 4]), htensor([4 4 4 4 4 4], 'TT'), [1 2 3 4], [1 5 6 2]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttt(htensor([4 4 4]), htensor([4 4 4 4]), [2 1 3], [1 3 4]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/ttv:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  ttv(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttv(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttv(x, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttv(x, {'a', 'b'});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttv(x, 1, 0.5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttv(x, 1, 8);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  ttv(x, 1, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/innerprod:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  innerprod(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  innerprod(x, x1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  innerprod(htensor([4 4]), htensor([3, 3]));
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  innerprod(x, {randn(4, 2), randn(5, 2)});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  innerprod(x, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  innerprod(x, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/cat:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  cat(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  cat(1, x, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  cat(1, x, x1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  cat(0.5, x, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  cat(8, x, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  cat(1, htensor([2 2 2]), htensor([2 2 3]));
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/nvecs:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  nvecs(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  nvecs(1, x, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  nvecs(x, 0.5, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  nvecs(x, 1, 0.5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/innerprod_mat:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  innerprod_mat(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  innerprod_mat(x, 1, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  innerprod_mat(x, x1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  innerprod_mat(x, x, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/mttkrp:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  mttkrp(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  mttkrp(1, 1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  mttkrp(x, 2, 3);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  mttkrp(x, {randn(3, 2), randn(4, 2), randn(5, 2), randn(6, 2)}, 0.5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  mttkrp(x, {randn(3, 2), randn(4, 2), randn(5, 2)}, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  mttkrp(x, {randn(3, 2), randn(4, 2), randn(5, 2), randn(5, 2)}, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  mttkrp(x, {randn(3, 2), randn(4, 2), randn(5, 2), randn(6, 3)}, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/power:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  power(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  power(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  power(x, 3);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/times:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  times(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x.*4;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  x.*x1;
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor([4 4 4]).*htensor([3 3 3]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/elem_mult:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  elem_mult(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  elem_mult(x, 4);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  elem_mult(x, x1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  elem_mult(htensor([4 4 4]), htensor([3 3 3]));
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  elem_mult(x, x, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/gramians_cp:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.gramians_cp(1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/gramians_sum:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.gramians_sum('a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.gramians_sum({'a', 'b'; 'c', 'd'});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.gramians_sum({'a', 'b'});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  htensor.gramians_sum({x, x1});
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/left_svd_gramian:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.left_svd_gramian('a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/left_svd_qr:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  htensor.left_svd_qr('a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/norm_diff:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  norm_diff(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  norm_diff(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  norm_diff(x, 'a');
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  norm_diff(x, [2 3 4]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  norm_diff(x, randn(size(x)), [2 3]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/full_block:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  full_block(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  full_block(x, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  full_block(x, [1 2; 1 2; 1 2; 1 100]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/full_mat:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  full_mat(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  full_mat(x, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  full_mat(x, 1, 2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  full_mat(x, [1 2 3 4], [1 2 3 4]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/apply_mat_to_vec:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  apply_mat_to_vec(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  apply_mat_to_vec(1, x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  apply_mat_to_vec(x, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  apply_mat_to_vec(x, x, 2);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  apply_mat_to_vec(htensor([3 3 3]), htensor([4 4 4]));
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(disp_message) fprintf('\n@htensor/apply_mat_to_mat:\n'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  apply_mat_to_mat(x);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  apply_mat_to_mat(1, x, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  apply_mat_to_mat(x, 1, 1);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  apply_mat_to_mat(x, x, 0.5);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  apply_mat_to_mat(htensor([4 4 4]), htensor([6 6 6]), [3 3 3]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end

try
  apply_mat_to_mat(htensor([6 6 6]), htensor([4 4 4]), [3 3 3]);
  disp('error not thrown');
catch err;
  if(disp_message) disp(err.message); end
end
