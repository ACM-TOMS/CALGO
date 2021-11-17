function B = horzcat(varargin)
%HORZCAT (overloaded)
%
% Author: Alexandre Felipe
% 2014, Dec, 8
%
%  Concatenates expressions side-by-side
B = varargin{1};
if(~isa(B, 'rolmipvar'))
  B = rolmipvar(B, '<>');
end
for i = 2:nargin
  h = homogenize(B, varargin{i});
  for j = 1:length(h{1}.data)
    B.data(j).exponent = h{1}.data(j).exponent;
    B.data(j).value = horzcat(h{1}.data(j).value, h{2}.data(j).value);
    B.opcode{j} = strcat(h{1}.opcode{j}, ',', h{2}.opcode{j});
  end
  B.vertices = h{1}.vertices;
  B.label = strcat(h{1}.label, ',', h{2}.label);
end
for j = 1:length(B.data)
  B.opcode{j} = strcat('[' , B.opcode{j}, ']');
end
B.label = strcat('[', B.label, ']');


