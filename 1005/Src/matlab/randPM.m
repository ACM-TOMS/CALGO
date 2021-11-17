% RANDPM  Uniform Park-Miller random numbers in [0,1)
%
%   This function uses the Park-Miller linear congruential generator described
%   in [1]. It enables programs written in C or Fortran produce exactly the same
%   sequence of random numbers as Matlab programs (to aid in testing and/or
%   debugging). Such a Fortran implementation is found in comments at the bottom
%   of this file.The built-in rand function is about 1000 times faster.

%
%   X = randPM() returns a uniformly distributed random number in [0,1)
%   X = randPM(N) returns an N by N matrix of uniform random numbers.
%   X = randPM(M,N) returns an M by N matrix of uniform random numbers.
%   X = randPM(N1,N2,...) returns an N1 by N2 by ... array.
%   randPM('init',S) sets the seed of the generator to S
%   randPM('init') sets the seed to 1234567
%
%   Reference
%   [1] S.K. Park and K.W. Miller (1988). "Randomn number generators: Good ones
%       are hard to find". Communications of the ACM 31 (10), 1192-1201.

function x = randPM(varargin)
  persistent a m q r fct seed
  if isempty(a)
    seed = 1234567;
    a = int32(16807);
    m = int32(2147483647);
    q = idivide(m,a);
    r = mod(m,a);
    fct = 1/(double(m)+1);
  end
  if nargin > 0 && strcmp(varargin{1}, 'init')
    if nargin < 2, seed = 1234567;
    else seed = varargin{2}; end
  elseif nargin > 0 && ischar(varargin{1})
    error('Illegal option');
  else
    if nargin == 0, varargin = {1,1}; end
    if nargin == 1, varargin{2} = varargin{1}; end
    n = prod(cell2mat(varargin));
    x = zeros(n,1);
    for i=1:n
      s1 = a*mod(seed,q);
      s2 = r*idivide(seed,q);
      if s1 > s2, seed = s1 - s2;
      else        seed = (m - s2) + s1; end
      x(i) = double(seed)*fct;
    end
    x = reshape(x, varargin{:});
  end
end

% real function rmd_srandpm(n)
%   integer n, seed, a, m, q, r, s1, s2
%   double precision fct
%   save fct, seed, a, m, q, r
% 
%   ! First run iff seed == 0 
%   if(seed .eq. 0) then 
%     seed = 1234567
%     a = 16807
%     m = 2147483647
%     q = m/a
%     r = mod(m,a)
%     fct = 1.0/(dble(m)+1.0)
%   endif
%   
%   ! If n > 0 then only set seed (and return)
%   if(n .gt. 0) then
%     seed = n 
%     return
%   endif
%   ! n <= 0, generate a pseudo-random number
%   
%   s1 = a*(mod(seed,q))
%   s2 = r*(seed/q)
%   
%   if(s1 .gt. s2) then
%     seed = s1-s2
%   else 
%     seed = (m - s2) + s1
%   endif
%   
%   rmd_srandpm = dble(seed)*fct
% 
% end function rmd_srandpm
