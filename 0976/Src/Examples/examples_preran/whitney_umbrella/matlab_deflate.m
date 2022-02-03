syms I Pi;
I=1i;
syms x y z;
X = [ x y z];
deg = [ 3];
f1=x^2-y^2*z;

F = [ f1];
J = jacobian(F,X);
R = nchoosek(1:1,1);
C = nchoosek(1:3,1);
r = size(R,1);
c = size(C,1);
count = 0;
OUT = fopen('deflation_polynomials','w');
for j = 1:r
  for k = 1:c
    A = simplify(det(J(R(j,:),C(k,:))));
    if A ~= 0
      count = count + 1;
		curr_eqn = char(2*A/prod(deg(R(j,:))));%
		i_locations = regexp(curr_eqn,'[\W\s]i[\W\s]');
		curr_eqn(i_locations+1) = 'I';
      fprintf(OUT, 'f_1_%d = %s;\n', count, curr_eqn);
    end;
  end;
end;
fclose(OUT);
OUT = fopen('deflation_polynomials_declaration','w');
fprintf(OUT, 'function ');
for j = 1:count
  fprintf(OUT, 'f_1_%d', j);
  if j == count    fprintf(OUT, ';\n');
  else    fprintf(OUT, ',');
  end;
end;

exit
