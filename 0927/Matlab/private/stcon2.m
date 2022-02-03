function [A_num2,B_num,C_num2] = stcon2()
      
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
      Rt21 = sqrt(21);
      A21 = 1/28 + 3*Rt21/1960;
      A22 = -Rt21/280;
      A23 = -32*Rt21/735;
      A24 = -23*Rt21/840;
      A25 = -1/28 + 3*Rt21/1960;

      A31 = 1/64;
      A32 = 7*Rt21/192;
      A34 = -7*Rt21/192;
      A35 = -1/64;

      A41 = 1/28 - 3*Rt21/1960;
      A42 = 23*Rt21/840;
      A43 = 32*Rt21/735;
      A44 = Rt21/280;
      A45 = -1/28 - 3*Rt21/1960;

      B1 = 1/20;
      B2 = 49/180;
      B3 = 16/45;

      C1 = 1/2 - Rt21/14;
      C2 = 1/2;
      C3 = 1/2 + Rt21/14;

      C12 = C1*C1;
      C22 = C2*C2;
      C32 = C3*C3;

      C16 = 6*(C1 - C12);
      C26 = 6*(C2- C22);
      C36 = 6*(C3 - C32);

      C123 = 3*C12 - 2*C1;
      C223 = 3*C22 - 2*C2;
      C323 = 3*C32 - 2*C3;

      C14 = 1 - 4*C1 + 3*C12;
      C24 = 1 - 4*C2 + 3*C22;
      C34 = 1 - 4*C3 + 3*C32;
     
      A_num2 =[A21 A22 A23 A24 A25
              A31 A32 0 A34 A35
              A41 A42 A43 A44 A45];
      
      B_num = [B1,B2,B3];
      
      C_num2 = [C1,C2,C3,C16,C26,C36,C123,C223,C323,C14,C24,C34];

end
