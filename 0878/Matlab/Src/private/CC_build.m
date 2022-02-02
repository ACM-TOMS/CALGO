%CC_BUILD  Calculate covariance of data and shocks
%
%  CC = CC_BUILD(A,C,n) finds the covariance between the vector x = [x1'...xn']'
%  and the vector [eps'...eps(n)']'. A is [A1...Ap] and C is {C0 C1...Cq}.
%  Used by varma_sim.

function CC = CC_build(A, C, n)
  r = size(C{1},1);
  q = length(C) - 1;
  p = size(A,2)/r;
  Ac = makecell(A);
  for j=q+1:p
    C{j+1} = zeros(r,r);
    for i=1:j, C{j+1} = C{j+1} + Ac{i}*C{j-i+1}; end
  end
  CC = cell(n,n);
  for i=1:n
    for j=1:n
      if 0<=i-j, CC{i,j} = C{i-j+1};
      else       CC{i,j} = zeros(r,r); end
    end
  end
  CC = cell2mat(CC);
end
