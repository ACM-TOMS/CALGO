function [def8,jc,problem] = df8cal_l(problem,h,nmsh,t,y,fty,ntol,ltol,tol,linear,A_num2,B_num,C_num2)

 %
%   Private function for twpbvpc
%
%   
%
%       Authors:
%
%       Jeff R. Cash 
%            (Department of Mathematics, Imperial College,  London, England.)
%       Davy  Hollevoet 
%            (Vakgroep Toegepaste Wiskunde en Informatica, Universiteit Gent, Belgium.)
%       Francesca Mazzia  
%            (Dipartimento di Matematica, Universita' di Bari, Italy)
%       Abdelhameed Nagy Abdo
%            (Dipartimento di Matematica, Universit\`a di Bari, Italy)
%            (Dept. of Mathematics, Faculty of Sciences, Benha  University,Egypt)
%            
%
    %  A_num2 =[A21 A22 A23 A24 A25
    %          A31 A32 0 A34 A35
    %          A41 A42 A43 A44 A45];
      
    %  B_num = [B1,B2,B3];
    %  
    %  C_num2 = [C1,C2,C3,C16,C26,C36,C123,C223,C323,C14,C24,C34];

index = 1:nmsh-1;
ncomp=problem.ncomp;
hb=repmat(h,ncomp,1);

tmp10 = zeros(3*ncomp,nmsh-1);

C16h = C_num2(4)./hb;
C26h = C_num2(5)./hb;
C36h = C_num2(6)./hb;

yn=y(:,1:end-1);
ynpo=y(:,2:end);
fn=fty(:,1:end-1);
fnpo=fty(:,2:end);



tmp4 = C16h.*(ynpo-yn)+ C_num2(7)*fnpo+ C_num2(10)*fn;
tmp5 = C26h.*(ynpo-yn)+ C_num2(8)*fnpo + C_num2(11)*fn;
tmp6 = C36h.*(ynpo-yn)+ C_num2(9)*fnpo + C_num2(12)*fn;

st1 = (yn+ynpo)/2;
st2 = A_num2(1,1)*fn + A_num2(1,5)*fnpo;
st3 = A_num2(2,1)*fn + A_num2(2,5)*fnpo;
st4 = A_num2(3,1)*fn + A_num2(3,5)*fnpo;

tn = t(1,1:end-1);
tc1 = tn + C_num2(1).*h;
tc2 = tn + C_num2(2).*h;
tc3 = tn + C_num2(3).*h;

dhold=zeros(3*ncomp,3*ncomp,nmsh-1);

for nit = 1:10
    tmp1 = st1(:,index) + hb(:,index).*(st2(:,index) + A_num2(1,2)*tmp4(:,index) + A_num2(1,3)*tmp5(:,index) + A_num2(1,4)*tmp6(:,index));
    tmp2 = st1(:,index) + hb(:,index).*(st3(:,index) + A_num2(2,2)*tmp4(:,index) + A_num2(2,4)*tmp6(:,index));
    tmp3 = st1(:,index) + hb(:,index).*(st4(:,index) + A_num2(3,2)*tmp4(:,index) + A_num2(3,3)*tmp5(:,index)+ A_num2(3,4)*tmp6(:,index));

    tmp7 = problem.f(tc1(index),tmp1)- tmp4(:,index);
    tmp8 = problem.f(tc2(index),tmp2)- tmp5(:,index);
    tmp9 = problem.f(tc3(index),tmp3)- tmp6(:,index);

    dfij=zeros(ncomp,ncomp,length(index));
    
    dfty=problem.df(tc1(index),tmp1);
   
    for i=1:length(index)
      dfij(:,:,i) = h(index(i))*dfty(:,:,i);  
    end
    dhold(1:ncomp,1:ncomp,index) = -A_num2(1,2)*dfij;
    dhold(1:ncomp,ncomp+1:2*ncomp,index) = -A_num2(1,3)*dfij;
    dhold(1:ncomp,2*ncomp+1:3*ncomp,index) = -A_num2(1,4)*dfij;

    
    dfty=problem.df(tc2(index),tmp2);

    for i=1:length(index)
      dfij(:,:,i) = h(index(i))*dfty(:,:,i);  
    end
    dhold(ncomp+1:2*ncomp,1:ncomp,index) = -A_num2(2,2)*dfij;
    dhold(ncomp+1:2*ncomp,ncomp+1:2*ncomp,index) = 0;
    dhold(ncomp+1:2*ncomp,2*ncomp+1:3*ncomp,index) = -A_num2(2,4)*dfij;
  
    
    dfty=problem.df(tc3(index),tmp3);

    for i=1:length(index)
      dfij(:,:,i) = h(index(i))*dfty(:,:,i);  
    end
    dhold(2*ncomp+1:3*ncomp,1:ncomp,index) = -A_num2(3,2)*dfij;
    dhold(2*ncomp+1:3*ncomp,ncomp+1:2*ncomp,index) = -A_num2(3,3)*dfij;
    dhold(2*ncomp+1:3*ncomp,2*ncomp+1:3*ncomp,index) = -A_num2(3,4)*dfij;

    
    for i=1:ncomp
        dhold(i,i,index)=dhold(i,i,index)+1;
        dhold(i+ncomp,i+ncomp,index)=dhold(i+ncomp,i+ncomp,index)+1;
        dhold(i+2*ncomp,i+2*ncomp,index)=dhold(i+2*ncomp,i+2*ncomp,index)+1;
    end

    for i=1:length(index)
        tmp10(:,index(i)) = dhold(:,:,index(i))\[tmp7(:,i);tmp8(:,i);tmp9(:,i)];
        tmp4(:,index(i)) = tmp4(:,index(i)) + tmp10(1:ncomp,index(i));
        tmp5(:,index(i)) = tmp5(:,index(i)) + tmp10(ncomp+1:2*ncomp,index(i));
        tmp6(:,index(i)) = tmp6(:,index(i)) + tmp10(2*ncomp+1:3*ncomp,index(i));
    end
    jc = 0;
    if (linear),  break, end
    indexz=zeros(1,nmsh-1);
    for i = 1: ntol
        ii = ltol(i);
        er = tol(i)./h;
        
        index = abs(tmp10(ii,:)) > er.*max(1,abs(tmp4(ii,:))) | abs(tmp10(ncomp+ii,:)) > er.*max(1,abs(tmp5(ii,:)))| abs(tmp10(2*ncomp+ii,:)) > er.*max(1,abs(tmp6(ii,:)));
        if any(index) % new modify
            jc = 1;
        end
        indexz(index)=1;
    end
    index=find(indexz == 1);
     
    if (jc == 0), break, end
end
if (linear || jc == 0)
   def8 = hb.*(B_num(1)*(fn + fnpo)+ B_num(2)*(tmp4 + tmp6)+ B_num(3)*tmp5)-ynpo+ yn;
else
  if (problem.debug), disp('no convergence of 8th order defcors'), end
  def8=0; 
end
if ~problem.vectorized 
     problem.NFUN = problem.NFUN + length(tc1(index))+ length(tc2(index)) + length(tc3(index));
   else
     problem.NFUN = problem.NFUN + 3;
end
end

