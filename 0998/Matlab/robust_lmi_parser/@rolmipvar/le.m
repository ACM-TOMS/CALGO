function LMIs = le(lhs, rhs)
%LE (overloaded less or equal)
%
% Author: Alexandre Felipe
% 2014, Dec, 8  LMIs = set([]);
  LMIs = [];
  h = homogenize(lhs, rhs);
  for i = 1:length(h{1}.data)
    LMIs = LMIs + (h{1}.data(i).value <= h{2}.data(i).value);
  end
end

