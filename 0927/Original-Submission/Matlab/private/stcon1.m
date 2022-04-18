function [A_num1,C_num1] = stcon1()
      
      %
%   Private function for twpbvpc
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
      Rt5 = sqrt(5);
      A21 = (6 + Rt5)/120;
      A22 = -Rt5/120;
      A23 = (-13 * Rt5)/120;
      A24 = (-6 + Rt5)/120;

      A31 = (6-Rt5)/120;
      A32 = (13 * Rt5)/120;
      A33 = Rt5 / 120;
      A34 = (-6 - Rt5)/120;

      C1 = (5 - Rt5)/10;
      C2 = (5 + Rt5)/10;

      C12 = C1*C1;
      C22 = C2*C2;

      C16 = 6*(C1 - C12);
      C26 = 6*(C2 - C22);

      C123 = 3*C12 - 2*C1;
      C223 = 3*C22 - 2*C2;

      C14 = 1 - 4*C1 + 3*C12;
      C24 = 1 - 4*C2 + 3*C22;
      
      A_num1 =[A21 A22 A23 A24
              A31 A32 A33 A34];
          
      C_num1 =[C1,C2,C16,C26,C123,C223,C14,C24];

 end
