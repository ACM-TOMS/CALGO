function  [def8,def6,jc,problem] = dfexcl_l(problem,h,nmsh,t,y,fty,ntol,ltol,tol,linear,A_num1,C_num1)
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
 %A_num1 =[A21 A22 A23 A24
 %         A31 A32 A33 A34];
          
 %  C_num1 =[C1,C2,C16,C26,C123,C223,C14,C24];
 
index = 1:nmsh-1;
ncomp=problem.ncomp;
hb=repmat(h,problem.ncomp,1);

tmp7 = zeros(2*ncomp,nmsh-1);

yn=y(:,1:end-1);          
ynpo=y(:,2:end);   

fn=fty(:,1:end-1);        
fnpo=fty(:,2:end);  

C16h = C_num1(3)./hb;
C26h = C_num1(4)./hb;

tmp3 = C16h.* (ynpo - yn) + C_num1(5)*fnpo + C_num1(7)*fn;
tmp4 = C26h.*(ynpo - yn) + C_num1(6)*fnpo + C_num1(8)*fn;
st1 = (yn+ynpo)/2;
st2 = A_num1(1,1)*fn + A_num1(1,4)*fnpo;
st3 = A_num1(2,1)*fn + A_num1(2,4)*fnpo;

tn = t(1,1:end-1);
tc1 = tn+ C_num1(1).*h;
tc2 = tn+ C_num1(2).*h;

dhold=zeros(2*ncomp,2*ncomp,nmsh-1);

for nit = 1:10 
    
    tmp1_f=st1 + hb.*(st2 + A_num1(1,2)*tmp3 + A_num1(1,3)*tmp4);
    tmp2_f=st1 + hb.*(st3 + A_num1(2,2)*tmp3 + A_num1(2,3)*tmp4);
    
    tmp1 = tmp1_f(:,index);
    tmp2 = tmp2_f(:,index);

    tmp5 = problem.f(tc1(index),tmp1)- tmp3(:,index);
    tmp6 = problem.f(tc2(index),tmp2)- tmp4(:,index);
    
    dfty=problem.df(tc1(index),tmp1);
    
    dfij=zeros(ncomp,ncomp,length(index));
    
    for i=1:length(index)
      dfij(:,:,i) = h(index(i))*dfty(:,:,i);
    end
    dhold(1:ncomp,1:ncomp,index) = -A_num1(1,2)*dfij;
    dhold(1:ncomp,ncomp+1:2*ncomp,index) = -A_num1(1,3)*dfij;
    
    
    dfty=problem.df(tc2(index),tmp2);
    
    for i=1:length(index)
      dfij(:,:,i) = h(index(i))*dfty(:,:,i); 
    end
    dhold(ncomp+1:2*ncomp,1:ncomp,index) = -A_num1(2,2)*dfij;
    dhold(ncomp+1:2*ncomp,ncomp+1:2*ncomp,index) =  -A_num1(2,3)*dfij;
    
     for i = 1:ncomp
           dhold(i,i,index) = dhold(i,i,index) + 1;
           dhold(i+ncomp,i+ncomp,index) = dhold(i+ncomp,i+ncomp,index) + 1;
     end
     
     for i=1:length(index)
       tmp7(:,index(i)) = dhold(:,:,index(i))\[tmp5(:,i);tmp6(:,i)];
       tmp3(:,index(i)) = tmp3(:,index(i)) + tmp7(1:ncomp,index(i));
       tmp4(:,index(i)) = tmp4(:,index(i)) + tmp7(ncomp+1:2*ncomp,index(i));
     end
     
     jc = 0;
     if (linear),  break, end

     indexz=zeros(1,nmsh-1);
     
     for i = 1: ntol
         ii = ltol(i);
         er = tol(i)./h; 
         
         indexz(abs(tmp7(ii,:)) > er.*max(1,abs(tmp3(ii,:)))  | abs(tmp7(ncomp+ii,:)) > er.*max(1,abs(tmp4(ii,:))))=1;
         
         if(any(indexz))
             jc=1;
         end
     end
     index=find(indexz == 1);
     
 if (jc == 0), break, end
end
if (linear || jc == 0)
    
    def6 = (hb./12).*(fn + 5*(tmp3 + tmp4)+fnpo) - ynpo + yn;
    for k=1:nmsh-1
        tmp7(:,k) = dhold(:,:,k)\[def6(:,k);def6(:,k)];
    end
    def8=tmp7(1:ncomp,:);
else  
    def6 = 0;
    def8 = 0;
   if (problem.debug),disp('no convergence of corrections'), end
end
if ~problem.vectorized 
     problem.NFUN = problem.NFUN + length(tc1(index))+ length(tc2(index));
   else
     problem.NFUN = problem.NFUN + 2;
end
end

     