function LMIs = gt(lhs, rhs)
%GT (overloaded greater than)
%
% Author: Alexandre Felipe
% 2014, Dec, 8  LMIs = set([])
  LMIs = [];
  h = homogenize(lhs, rhs);
%   if(isa(rhs, 'rolmipvar'))
    for i = 1:length(h{1}.data)
      LMIs = LMIs + (h{1}.data(i).value > h{2}.data(i).value);
    end
%   end
end
