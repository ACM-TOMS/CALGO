function [y, bkerr, pjcnd, job ] = ...
              pejroot( f, y0, l, noi, tol)
% PEJROOT calculates (multiple) roots of polynomial 
%
% This code is freely released for research exchange only. The author 
% is not responsible for any demage caused by using this code.
%
% Calling syntax: The simplest way to call is
%  
%   © z = pejroot(f,y,m)
%
% where the input items are
%          f --- (row vector) the polynomial
%          y --- (row vector) the initial iterate (of the roots)
%          m --- (row vector) the multiplicity structure
%
% For more advanced usage:
%
%	© [z, e, c] = pejroot(f, y, l, noi, tol)
%
% The output 
%          z --- (matrix) distinct roots (1st column) and corresponding
%                         multiplicities (2nd column)
%          e --- (scalar) the backward error of the roots
%          c --- (scalar) the pejorative condition number
%
%
%  Additional input parameters
%  
%        noi --- (integer)    the number of iteration allowed
%                               (default: 10)
%        tol --- (numeric)    the error tolerance
%                               (default: 1.0d-8)

   if nargin == 3, noi = 10; tol = 1.e-6; end;
   
   m = length(l);    % number of variables
   n = sum(l);       % number of equations
   
   if abs(f(1)) ~= 0     % make the polynomial monic
       f = f / f(1);
   else, jj = find(f);  j = min(jj); f = f(j:length(f));
   end;
   
   if length(f) ~= n+1   % exit if the input is wrong
       fprintf('In put error');  job = 1; return;
   end;
       
   % sort initial values to enhance accuracy
   % it is interesting, sorting really improves accuracy
   [yy,jj] = sort(1./abs(y0));   y0 = y0(jj); l = l(jj);

   job = 0;               % initialze job
   y = conj(y0');         % pass the initial iterate
   h = conj(f(2:n+1)');   % make the RHS
   delta = zeros(1,noi);  % space for sizes of the correction
   bcker = zeros(1,noi);  % space for backward errors

   w = ones(n,1);         % set weight
   for j = 1:n
       if abs(f(j+1)) >= 1.0 
           w(j) = abs(1.0/f(j+1));
       end
   end;

   for k = 1:noi

      % evaluate the coefficient operator and its Jacobian
      Df = zeros(n,m);       % open the space for A
      
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
         g = conv(-l(j)*g,ff);  Df(:,j) = conj(g');
      end
      g = conv(-g/l(m), [1,-y(m)]);   % evalulate the coef. operator
      c = conj(g(2:n+1)');
      
      b = c - h;   % RHS of the polyn. system
 
      d = scalsq(Df,b,w);
      delta(k) = norm(d,2);   bcker(k) = norm(w.*b,inf);

      if delta(k) < tol, job = 1; end     % convergence criterion 1
      if k > 1 
         if delta(k) > delta(k-1) & bcker(k) > bcker(k-1)
            if job == 1, bkerr = bcker(k-1);   break; end;
         elseif delta(k) < delta(k-1)                 % criterion 2
            if delta(k)^2/(delta(k-1)-delta(k)) < tol, 
                job = 1; 
            end
         end
      end
      
      y = y - d;    % correct the roots
  
      bkerr = bcker(k); 
      A = Df;
   end

   if job == 1
      s = svd(A); pjcnd = 1/s(m);               % get pej. cond. number
      [ll,jj] = sort(l);  y = [y(jj),l(jj)'];   % sort by multiplicities
   else
      pjcnd = 0;
   end;
