function P=SortAndRemoveDoubles(P);
% SORTANDREMOVEDOUBLES sort points and remove duplicates
%   P=SortAndRemoveDoubles(P); orders polygon corners in P counter
%   clock wise and removes duplicates

ep=10*eps;                           % tolerance for identical nodes
m=size(P,2);
if m>0
  c=sum(P')'/m;                      % order polygon corners counter
  for i=1:m                          % clockwise
    d=P(:,i)-c; ao(i)=angle(d(1)+sqrt(-1)*d(2));
  end;
  [tmp,id]=sort(ao);
  P=P(:,id);
  i=1;j=2;                           % remove duplicates
  while j<=m
    if norm(P(:,i)-P(:,j))>ep
      i=i+1;P(:,i)=P(:,j);j=j+1;
    else
      j=j+1;
    end;
  end;
  P=P(:,1:i);
end;
