function [p,z] = toh06a
%
%  generalization of K.C. Toh and L. N. Trefethen
% 
   p = ones(1,21);
     z = [  0.95557280578614 + 0.29475517441090i
  0.95557280578614 - 0.29475517441090i
  0.82623877431600 + 0.56332005806362i
  0.82623877431600 - 0.56332005806362i
  0.62348980185873 + 0.78183148246803i
  0.62348980185873 - 0.78183148246803i
  0.36534102436640 + 0.93087374864420i
  0.36534102436640 - 0.93087374864420i
  0.07473009358642 + 0.99720379718118i
  0.07473009358642 - 0.99720379718118i
 -0.98883082622513 + 0.14904226617617i
 -0.98883082622513 - 0.14904226617617i
 -0.90096886790242 + 0.43388373911756i
 -0.90096886790242 - 0.43388373911756i
 -0.73305187182983 + 0.68017273777092i
 -0.73305187182983 - 0.68017273777092i
 -0.50000000000000 + 0.86602540378444i
 -0.50000000000000 - 0.86602540378444i
 -0.22252093395631 + 0.97492791218182i
 -0.22252093395631 - 0.97492791218182i];
    z = [z,ones(20,1)];
   if norm(imag(z(:,1))) == 0 
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.15f \t \t \t %3g \n', z');
   else
        fprintf('                 roots ')
        fprintf('   \t\t\t\t\t\t     multiplicities\n');
        fprintf('\n');
        fprintf('%22.15f + %22.15f i \t \t %3g \n', ...
            [real(z(:,1)),imag(z(:,1)),z(:,2)]');

   end;       
