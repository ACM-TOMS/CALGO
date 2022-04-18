function disp(f)

% $Id: disp.m,v 1.3 2001/07/20 14:03:14 driscoll Exp $ 

s = char(f);
if isstr(s)
  disp(s)
elseif iscell(s)
  fprintf('\n  SC %s:\n\n',class(f));
  for n = 1:length(s)
    disp(s{n})
  end
end
