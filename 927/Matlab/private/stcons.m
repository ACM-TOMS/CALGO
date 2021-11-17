function [c_alp, c_bet,c_abc]=stcons()
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
     
 
%  stcons computes constants needed in integration formulae
%  
%     c_abc = [ a1 b1 c1 d1 0  0 0  0  0
%               a2 b2 c2 d2 e2 f2 0  0  0 
%               a3 b3 c3 d3 e3 f3 p3 q3 0 
%               a4 0  c4 0  e4 0  p4 0  x4
%               a5 b5 c5 d5 e5 f5 0  0  0
%               a6 b6 c6  0  0  0 0  0  0] ;  
%
      one = 1.0d+0;
      four = 4.0d+0;
      two = 2.0d+0 ;
      five = 5.0d+0;
      three = 3.0d+0;
      half = 0.5d+0;
      fourth = 0.25d+0;   
     
           

      alp1 = one/8.0d+0;
      alp1sq = one/64.0d+0;
      alp2sq = five/28.0d+0;
      alp3sq = five/84.0d+0;
      alp2 = sqrt(alp2sq);
      alp3 = sqrt(alp3sq);
      bet0 = one/6.0d+0 - (10.0d+0 - 28.0d+0*(alp2sq + alp3sq))/...
              (105.0d+0*(one - four*alp2sq)*...
              (one - four*alp3sq));
      bet2 = -(28.0d+0*alp3sq - three)/...
              (1680.0d+0*alp2sq*...
              (one - four*alp2sq)*(alp2sq - alp3sq)); 
      bet3 = (28.0d+0*alp2sq - three)/...
               (1680.0d+0*alp3sq*...
               (one - four*alp3sq)*(alp2sq - alp3sq));
      bet4 = half - bet0 - bet2-bet3;
      
      a1 = half*(one - alp1)*(two*alp1 + one)*...
               (two*alp1 + one) ;
      b1 = half*(one + alp1)*...
               (two*alp1 - one)*(two*alp1 - one); 
      c1 = half*(alp1sq - fourth)*(two*alp1 + one) ;
      d1 = half*(alp1sq - fourth)*(two*alp1 - one); 
      uu = alp2*((four*alp2sq - one)^2)/...
            ((four*alp1sq - one)*(20.0d+0*alp1*alp1 - one)) ;
      vv = ((four*alp2sq - one)^2)/...
           (16.0d+0*alp1*(four*alp1sq - one)) ;
      e2 = half*(uu + vv) ;
      f2 = half*(uu - vv) ;
      rr = half*(alp2*(four*alp2sq - one) +...
           (one - 12.0d+0*alp1sq)*(e2 + f2));
      ss = fourth*(four*alp2sq - one) - two*alp1*(e2 - f2); 
      c2 = half*(rr + ss) ;
      d2 = half*(rr - ss) ;
      ww = two*(alp2 - (c2 + d2 + e2 + f2));
      b2 = half*(one - ww) ;
      a2 = half*(one + ww) ;
      z1 = (three - 28.0d+0*alp3sq)/...
              (1680.0d+0*alp2*(four*alp2sq - one)*...
              (alp2*alp2 - alp3sq)*bet3) ;
      z2 = one/(105.0d+0*alp3*bet3*...
         (20.0d+0*alp2sq - one)*(four*alp2sq - one)) ;
      p3 = half*(z1 + z2); 
      q3 = half*(z2 - z1) ;
      u1 = (alp3*((four*alp3sq - one)^2)-...
            (p3 + q3)*(20.0d+0*alp2sq - one)...
            *(four*alp2sq - one))/...
             ((four*alp1sq - one)*(20.0d+0*alp1sq - one)) ;
      v1 = (alp3sq*(one - two*alp3sq)-...
              two*alp2*(one - four*alp2sq)*(p3 - q3)...
             -one/8.0d+0)/(two*alp1*(one - four*alp1sq));
      e3 = half*(u1 + v1) ;
      f3 = half*(u1 - v1) ;
      r1 = half*(alp3*(four*alp3sq - one) +...
              (e3 + f3)*(one - 12.0d+0*alp1sq) + ...
              (p3 + q3)*(one - 12.0d+0*alp2sq));
      s1 = alp3sq - fourth - two*alp1*(e3 - f3)...
              - two*alp2*(p3 - q3); 
      c3 = half*(r1 + s1); 
      d3 = half*(r1 - s1); 
      w1 = two*(alp3 - (c3 + d3 + e3 + f3 + p3 + q3));
      a3 = half*(one + w1); 
      b3 = half*(one - w1) ;
      a4 = half;
      p4 = 0.0d+0 ;
      x4 = (three - 28.0d+0*alp2sq)/...
             (3360.0d+0*alp3*bet4*(four*alp3sq - one) ...
             *(alp3sq - alp2sq));
      e4 = (0.125d+0 + four*alp2*p4*...
              (one - four*alp2sq) +...
              four*alp3*x4*(one - four*alp3sq))/...
              (four*alp1*(four*alp1sq - one)) ;
      c4 = -(0.125d+0 + two*alp1*e4 + two*alp2*p4 + ...
              two*alp3*x4);
      a5 = five/32.0d+0 ;
      b5 = 27.0d+0/32.0d+0;
      c5 = 9.0d+0/64.0d+0 ;
      d5 = three/64.0d+0 ;
      e5 = five/24.0d+0; 
      f5 = two/three;
      a6 = 7.0d+0/90.0d+0 ;
      b6 = 16.0d+0/45.0d+0;
      c6 = two/15.0d+0; 
    
      c_alp = [alp1,alp2,alp3];
      c_bet = [bet0,bet2,bet3,bet4];
      c_abc = [ a1 b1 c1 d1 0  0 0  0  0
               a2 b2 c2 d2 e2 f2 0  0  0 
               a3 b3 c3 d3 e3 f3 p3 q3 0 
               a4 0  c4 0  e4 0  p4 0  x4
               a5 b5 c5 d5 e5 f5 0  0  0
               a6 b6 c6  0  0  0 0  0  0] ;  
end

  