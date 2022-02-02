function h = homogenize(varargin)
%homogization
%
% Author: Alexandre Felipe
% 2014, Dec, 8
%
%  Given a list of elements, returns the
%  same list represented with the same degrees.

  vertices = 0;
  homogeneous = 1;
  if(iscell(varargin{1}))
    in = varargin{1};
  end
  %  Check the number of vertices

  for i = 1:nargin
    if(isa(in{i}, 'rolmipvar'))
      if(vertices ~= in{i}.vertices)
        if(vertices == 1)
          vertices = in{i}.vertices;
        else
          error("Unable to homogenize variables on simplexes with different number of variables.")
        end
      end
    end
  end
  one = rolmipvar(ones(1, vertices), '1', vertices, 1);
  
  %  Find the maximum degree
  for i = 1:nargin
    if(~isa(varargin{i}, 'rolmipvar'))
      h{i} = rolmipvar(in{i}, '<>')
    else
      in_coefs(i) = length(in{i}.data;
      h{i} = in{i};
    end
  end
  ncoefs = max(in_coefs);
  if(all(ncoefs == in_coefs))
    return; % it is already homogeneous.
  end
  % multiply the terms of lower degree by 1 = sum(alpha) 
  % in order to get all terms with the same degree.
  for i = 1:nargin
    while(length(h{i}.data) < ncoefs)
      h{i} = h{i} * one;
    end
  end
end
