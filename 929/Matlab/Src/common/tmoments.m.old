function TMom = tmoments (Mom,i);
%TMOMENTS  Moments of translated scaling functiones
%
%  TMom = tmoments (Mom,l);
% 
%  Input
%    Mom: Vector with moments obtained from moment.m
%    l:   Translation parameter
%
%  Output
%    TMom: The translated moment
%
%
%  See also moment, conn.


M0 = 1;
pmax = length(Mom);
TMom = zeros (size(Mom));

for j = 1:pmax
  for k = 0:j
    if k==0 
      Mk = M0;
    else
      Mk = Mom(k);
    end;  
    
    bin = nchoosek(j,k);
    TMom(j) = TMom(j) + bin * Mk * i^(j-k);
  end
end

