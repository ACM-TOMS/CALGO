% MAKECELL  Change to cell array
%
%   C = MAKECELL(C) changes C from [C1 C2...Cm] to {C1 C2...Cm} where each Ci is
%   square. If C is empty the empty cell array, {}, is returned.

function Cout = makecell(C)
  if isempty(C), Cout = {}; return, end
  r = size(C,1);
  m = size(C,2)/r;
  Cout = cell(1, m);
  for i=1:m
    Cout{i} = C(:, (i-1)*r+1:i*r);
  end
end
