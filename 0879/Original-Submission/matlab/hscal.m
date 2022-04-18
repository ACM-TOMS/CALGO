function v = hscal(u)
%
%  v = hscal(u)
%
%  scales the nonzero vector u so that its norm is sqrt(2).
%

   nrm = norm(u);
   if nrm == 0
      error('Error in hnorm: zero argument.');
   end
   v = (sqrt(2)/nrm)*u;
return
