function Mom = moments (hk, pmax);
%MOMENTS  Compute moments of scaling function 
%
%  Mom = moments (hk, pmax); 
%
%  Input:
%    hk:   Filter coefficients (e.g. obtained from wfilters)
%    pmax: Number of moments to be computed
%
%  Output:
%    Mom:  Vector of moments from 1 to pmax
%
% Use transmoments.m to obtain moments of translated scaling function
%
% The zeroth moment is always 1 and is therefore not stored
%
% See also transmoments, conn.
N = length (hk);
M0  = 1;
Mom = zeros (1, pmax);
for p = 1:pmax
  fact = sqrt(2)/2/(2^p-1);  
  for k = 0:p-1
    if k==0 
      Mk = M0;
    else
      Mk = Mom(k);
    end;  
    dM = 0;              %"Discrete moment"
    for i=0:N-1
    dM = dM + hk(i+1) * i^(p-k);
    end;
bin    = nchoosek(p,k);
 Mom(p) =  Mom(p)+bin * Mk * dM; % dbm
  end
  Mom(p)=fact*Mom(p);  
end




