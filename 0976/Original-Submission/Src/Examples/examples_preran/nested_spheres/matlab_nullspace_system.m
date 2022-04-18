syms I Pi;
syms x y z u v;
my_var_names = [ x y z u v];

syms pi_1_2 pi_1_3 pi_1_4 pi_1_5 pi_1_6 pi_2_2 pi_2_3 pi_2_4 pi_2_5 pi_2_6
proj = [[ pi_1_2 pi_1_3 pi_1_4 pi_1_5 pi_1_6 ]; [ pi_2_2 pi_2_3 pi_2_4 pi_2_5 pi_2_6 ]; ];

num_projections = 2;
num_jac_equations = 5;
num_randomized_eqns = 3;
target_crit_codim = -1;
r_squared=x^2+y^2+z^2;
f1=r_squared-1-u^2;
f2=r_squared-4+v^2;
crit=-4*u*v;

syms r_1_1 r_1_2 r_1_3 r_2_1 r_2_2 r_2_3 r_3_1 r_3_2 r_3_3 ;

F_orig = [ f1; f2; crit;]; %collect the functions into a single matrix
F_rand = F_orig; % no need to randomize
J = [jacobian(F_rand,my_var_names); proj];  %compute the transpose of the jacobian
                                           %concatenate the projections

new_eqn = det(J); 

OUT = fopen('derivative_polynomials_declaration','w'); %open the file to write
fprintf(OUT, 'function ');
for ii=1:num_randomized_eqns
  fprintf(OUT, 'f%i',ii);
  if ii~=num_randomized_eqns
    fprintf(OUT,', ');
  else
    fprintf(OUT,';\n');
  end %re: if
end

for ii=1:num_randomized_eqns
  fprintf(OUT,'f%i = %s;\n',ii,char(F_rand(ii)));
end

fprintf(OUT, 'function der_func;\n');
fprintf(OUT,'der_func = %s;\n',char(new_eqn));

exit %quit the script
