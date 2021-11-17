%    This is a matlab-octave script file which reduces
%    a general square real matrix to small band Hessenberg form.
%    It requires as input
%
%	a   --  square matrix.
%       n   --  square matrix a is n by n.
%      tol  -- a nonnegative number which is the bound
%                      on the multipliers.  
%
%    output b -- banded Hessenberg matrix
%    output z -- accumulated similarity transformation z 
%
%    This routine is not very efficient and has many more
%    variables than required. It is designed to  check
%    correct execution of the FORTRAN algorithm BHESS
%    and auxiliary algorithms 
%
%    This routine explicitly computes the accumulated similarity
%    transformation z
%
%            b = inv(z) * a * z 
%
%     allowing verifcation of BHESS and also auxiliary routines 
%     packaged as auxiliary FORTRAN routines.
% 
%     To use the same matrix in both the Fortran and matlab
%     version, run the matlab version first, then use
%
%     save BHESS.IN tol n a
%
%     Edit BHESS.IN to get rid of nonnumeric lines.  
%
%     Then run the BHESS driver FORTRAN executable
%     and take the option of reading the input matrix
%     from BHESS.IN
%  
%     You can then compare the matrices returned and the
%     results of multiplying a vector by z or inv(z).  
%  
%     Results agree to more digits if 
%
%    format long
% 
%     is used inside matlab or octave before creating the BHESS.In
%     file.     
%
%
b=a; 
ltot=eye(n);
j=1
sn=eye(n);
sv=eye(n);
sb=eye(n);
ao=ones(n);
kr=ao(1:n,1);      
%  loop over columns , zeroing each in turn below subdiagonal 
for i=2:n-1,
    l=eye(n);
    xb=b(i:n,i-1);
    kt=1;
    kc=j;
%  Test which if any row can be eliminated. 
    while( (kc<i)&(kt==1) )
      if(kr(kc)==1)
        xrr=b(kc,i:n);
        cs=xrr*xb/norm(xrr)/norm(xb);
        kc=kc;
        if (max(abs(cs))>(tol^(-1)/(n-i+1)^1.0))
% This row can be eliminated, so find the right pivot.  
          inprod = xrr*xb;
          [yb,k1]=max(abs(xrr));
          pivr = k1;
          xrr2=xrr;
          xrr2(k1)=0; 
          [yb2,k2]=max(abs(xrr));
          [ybc,k3]=max(abs(xb));
          xb2=xb; 
          xb2(k3)=0; 
          pivc=k3; 
          [ybc2,k4]=max(abs(xb2));
          pivp=1; 
          minmlt = 1.e10;
          for ii=1:n-i+1,
            v1=xrr(ii);
            if(v1 == 0.0)
              maxnr = 1.e10;
            elseif (ii == pivr)
              maxnr = abs(yb2/v1);
            else maxnr=abs(yb/v1);
            end
            if (ii == pivc)
              maxnc=abs(ybc2*v1/inprod);
            else 
              maxnc=abs(ybc*v1/inprod) ;
            end 
            mxdiag = abs(v1*xb(ii)/inprod);
            temp = maxnr;
            if (temp < maxnc) temp = maxnc; end;
            if (temp < mxdiag) temp = mxdiag; end; 
            if (temp < minmlt) 
              minmlt = temp ;
              pivp = ii;
            end
          end 
          k = pivp; 
          kr(kc) = 0; 
          kt = 0; 
        end
      end
      kc=kc+1;
    end
    if(kt==1)
      [yb,k]=max(abs(xb)); % this is the case that no row can be
%			eliminated so choose k to swap
% 			largest column element to subdiagonal.  
    end
    p=eye(n);
    p(i,i)=0;
    p(i+k-1,i+k-1)=0;
    p(i,i+k-1)=1;
    p(k+i-1,i)=1;
% permute ith and kth rows and columns 
    b=p*b*p;
    ltemp=l;
    ltempn=l*ltemp;
%  A row can be eliminated with the pivot chosen above 
    if (kt==0)
      xr=b(kc-1,i+1:n)/b(kc-1,i);
      r=eye(n);
      r(i,i+1:n)=-b(kc-1,i+1:n)/b(kc-1,i);
% Spending O(n^3) flops for one row elimination; 
%   is not very efficient.  
      b=inv(r)*b*r;
      ltempn=inv(r)*ltempn;
      ltemp =inv(r)*ltemp;
    end
    l(i+1:n,i)=-b(i+1:n,i-1)/b(i,i-1);
% and O(n^3) flops for one column elimination.  
%   Also inefficient. 
    b=l*b*inv(l);
% those who want a really detailed snap shot of the
% process may enjoy getting the svd of the 
% accumulated transformation on each step.  
%  This is of course even more ineffficient.  
    if (kt==0)
      ltot=l*inv(r)*p*ltot;
      sn(:,i)=svd(ltempn);
    else
      ltot=l*p*ltot;
      sn(:,i)=svd(ltempn);
    end
    sv(:,i)=svd(ltot);
    sb(i:n,i)=svd(b(i:n,j:n));
    x=sn(:,i-1:i);y=sv(:,i-1:i);
    z=sb(:,i)
    pause
%   pr=prod(sv)
    co=sv(1,i)/sv(n,i);
end;
b=b
z=inv(ltot)
size=norm(b)
size2=norm(ltot)
co=cond(ltot)  % you may want to compare this to the estimated
% condition number in option 4 of the driver routine.  
size3=norm(a/ltot)
s=0
for i=1:n,
  for j=1:n,
    s=max(s,abs(ltot(i,j)));
  end;
end

