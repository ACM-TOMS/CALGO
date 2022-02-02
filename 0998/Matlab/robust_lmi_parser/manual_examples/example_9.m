clear;

A1 = eye(2);
A2 = 2*eye(2);
A = [A1 A2];
poly_A = rolmipvar(A,'A',2,1);

B1 = 3*eye(2);
B2 = 4*eye(2);
B = [B1 B2];
poly_B = rolmipvar(B,'B',2,1);

fprintf(1,'A+B:\n');
resul = poly_A + poly_B

fprintf(1,'\nA*B:\n');
resul = poly_A*poly_B

fprintf(1,'\nA+C:\n');
C = 5*eye(2);
resul = poly_A + C


fprintf(1,'\nA*C:\n');
resul = poly_A*C