function cond = ...
                 spcond(z)
% SPCond calculates the struture preserving condition number
% of roots z(:,1) with multiplicity structure z(:,2)
% 
%   Input:  z --- mx2 matrix, z(:,1) is the root vector
%                             z(:,2) is the multiplicity structure
%
%  Output:  cond -- the structure preserving condition number
%
% Syntax:  >> cond = spcond(z)
%
   m = size(z,1);   % number of distinct roots
   n = sum(z(:,2)); % degree of the polynomial
   
   % sort
   [yy,jj] = sort(1./abs(z(:,1)));   y = z(jj,1); l = z(jj,2)';
   
   ff = [1];  ll = l - ones(1,m);    % form the base polynomial
   mm = max(ll); lll = ll;
   for kk = 1:mm
       id = find(lll>0);
       ff = conv(ff,poly(y(id)));
       lll = lll - 1;
   end;
   
   for j = 1:m               % remaining polynomial multiplication
       if m > 1
           id = [1:(j-1),(j+1):m]; g = poly(y(id));
       else
            g = [1];
       end;
       g = conv(-l(j)*g,ff);  J(:,j) = conj(g');
   end
   g = conv(-g/l(m), [1,-y(m)]);  % evalulate the coef. operator
   c = conj(g(2:n+1)');
  
   for j = 1:n                    % scale J
       if abs(c(j)) >= 1.0 
           s = abs(1.0/c(j));
           J(j,:) = J(j,:)*s;
       end
   end;
   
   [s,x] = zminsv(J,1.0e-16);       % finding smallest singular value
   
   cond = 1/s;
