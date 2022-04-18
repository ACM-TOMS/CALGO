function p=RemoveDuplicates(p)
% REMOVEDUPLICATES removes multiple entries
%   p=RemoveDuplicates(p) removes duplicate entries from the vector
%   p.

p=sort(p);
i=1;j=2;                          % remove duplicates
while j<=length(p)
  if p(i)<p(j)
    i=i+1;p(i)=p(j);j=j+1;
  else
    j=j+1;
  end;
end;
p=p(1:i);
