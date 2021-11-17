function [def8,problem]=df8cal_m(problem,h,t,y,fty,c_alp,c_bet,c_abc)
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
    fc1=0.625;
    fc2=0.375;
     
    hb=repmat(h,problem.ncomp,1);
    % todo: pass hb around?
    
    yn=y(:,1:end-1);
    ynpo=y(:,2:end);
    fn=fty(:,1:end-1);
    fnpo=fty(:,2:end);
    
    
    tn=t(1,1:end-1);
    tnpof=tn+fc1*h;
    tnpoh=tn+fc2*h;
   
    
    tmp1 = c_abc(1,1)*ynpo+c_abc(1,2)*yn+hb.*(c_abc(1,3)*fnpo+c_abc(1,4)*fn);
    tmp2 = c_abc(1,2)*ynpo+c_abc(1,1)*yn-hb.*(c_abc(1,3)*fn+c_abc(1,4)*fnpo);
    
    tmp3=problem.f(tnpof,tmp1);
    tmp4=problem.f(tnpoh,tmp2);
     
    tmp1 = c_abc(2,1)*ynpo+c_abc(2,2)*yn+hb.*(c_abc(2,3)*fnpo+c_abc(2,4)*fn+ ...
           c_abc(2,5)*tmp3+c_abc(2,6)*tmp4);
    tmp2 = c_abc(2,2)*ynpo+c_abc(2,1)*yn-hb.*(c_abc(2,3)*fn+c_abc(2,4)*fnpo+ ...
           c_abc(2,5)*tmp4+c_abc(2,6)*tmp3);

       
     tnptf=tn + (0.5 + c_alp(2))*h;
     tnpth=tn + (0.5 - c_alp(2))*h;
     
     tmp5=problem.f(tnptf,tmp1);
     tmp6=problem.f(tnpth,tmp2);
  
     tmp1 = c_abc(3,1)*ynpo+c_abc(3,2)*yn+hb.*(c_abc(3,3)*fnpo+c_abc(3,4)*fn+ ...
            c_abc(3,5)*tmp3+c_abc(3,6)*tmp4 + ...
            c_abc(3,7)*tmp5+c_abc(3,8)*tmp6) ;
     tmp2 = c_abc(3,2)*ynpo+c_abc(3,1)*yn-hb.*(c_abc(3,3)*fn+c_abc(3,4)*fnpo+ ...
           c_abc(3,5)*tmp4+c_abc(3,6)*tmp3+ ...
           c_abc(3,7)*tmp6+c_abc(3,8)*tmp5) ;
       
     tnptf=tn + (0.5 + c_alp(3))*h;
     tnpth=tn + (0.5 - c_alp(3))*h;
     
     tmp7=problem.f(tnptf,tmp1);
     tmp8=problem.f(tnpth,tmp2);

     tmp1 = c_abc(4,1)*(ynpo+yn)+hb.*(c_abc(4,3)*(fnpo-fn)+ ...
            c_abc(4,5)*(tmp3-tmp4) + ...
            c_abc(4,9)*(tmp7-tmp8)) ;
     
        
      tnpth=tn + 0.5*h;
      tmp2=problem.f(tnpth,tmp1);
     
        
      
      def8 = hb.*(c_bet(1)*(fn+fnpo)+c_bet(2)*(tmp5+tmp6)+c_bet(3)*(tmp7+tmp8)+...
                  2*c_bet(4)*tmp2)-ynpo+yn;
              
    if ~problem.vectorized 
        problem.NFUN = problem.NFUN + length(tnpof) +  length(tnpoh) + length(tnptf) + length(tnpth) ;
    else
        problem.NFUN = problem.NFUN + 4;
    end
   
end