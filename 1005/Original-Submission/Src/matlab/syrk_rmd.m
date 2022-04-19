% SYRK_RMD  Reverse mode derivative of Blas syrk operation

% Note: Had to rename sym to symm due to bug in Matlab R2018a
function Aa = syrk_rmd(A, Aa, La, trans)
  trsp = nargin > 3 && strcmpi(trans, 't');
  if trsp
    Aa = Aa + A*(symm(La) + dg(La));  
  else
    Aa = Aa + (symm(La) + dg(La))*A;   
  end
end
