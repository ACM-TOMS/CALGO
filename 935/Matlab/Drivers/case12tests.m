
addpath('..')

rho_vec = [0.001 0.01 0.1 1 10 100 1000];
tau_vec = [0.0011 0.011 0.11 1.1 11 101 1001];

abserr = 1e-10;relerr = 1.0e-11;
u = 1; a = 0; b = 0;type='JJ';
fx = @(x)exp(-u*x);
for i=1:length(rho_vec)
    for j=1:length(tau_vec)
        [results(i,j),rel_err(i,j), neval(i,j)]=IIPBF(fx,rho_vec(i),tau_vec(j),a,b,abserr,relerr,type);
        temp1=1/(sqrt(u^2 + (rho_vec(i)+tau_vec(j))^2));
        temp=2/pi*temp1;
        answer= temp*ellipk(2*sqrt(rho_vec(i)*tau_vec(j))*temp1);
        act_err(i,j) = (abs(results(i,j) - answer));
        act_err_1(i,j) = act_err(i,j)/(10.^floor(log10(act_err(i,j))));
        act_err_2(i,j) = floor(log10(act_err(i,j)));
        if (act_err(i,j)==0)
            act_err_1(i,j) = 0;
            act_err_2(i,j) = 0;
        end
    end
end

fprintf('\n\\multicolumn{9}{c}{Case 12} \\\\\n');
fprintf('       &        &               \\multicolumn{7}{c}{$\\rho$} \\\\\n');
fprintf('       &        & %s & %s & %s & %s & %s & %s & %s \\\\ \\cline{3-9}\n', num2str(rho_vec(1)), num2str(rho_vec(2)), num2str(rho_vec(3)), num2str(rho_vec(4)), num2str(rho_vec(5)), num2str(rho_vec(6)), num2str(rho_vec(7)));
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(1)), act_err_1(1,1), ['x$10^{' num2str(act_err_2(1,1)) '}$' ],  act_err_1(2,1), ['x$10^{' num2str(act_err_2(2,1)) '}$' ], act_err_1(3,1), ['x$10^{' num2str(act_err_2(3,1)) '}$' ], act_err_1(4,1), ['x$10^{' num2str(act_err_2(4,1)) '}$' ], act_err_1(5,1), ['x$10^{' num2str(act_err_2(5,1)) '}$' ], act_err_1(6,1), ['x$10^{' num2str(act_err_2(6,1)) '}$' ], act_err_1(7,1), ['x$10^{' num2str(act_err_2(7,1)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(2)), act_err_1(1,2), ['x$10^{' num2str(act_err_2(1,2)) '}$' ],  act_err_1(2,2), ['x$10^{' num2str(act_err_2(2,2)) '}$' ], act_err_1(3,2), ['x$10^{' num2str(act_err_2(3,2)) '}$' ], act_err_1(4,2), ['x$10^{' num2str(act_err_2(4,2)) '}$' ], act_err_1(5,2), ['x$10^{' num2str(act_err_2(5,2)) '}$' ], act_err_1(6,2), ['x$10^{' num2str(act_err_2(6,2)) '}$' ], act_err_1(7,2), ['x$10^{' num2str(act_err_2(7,2)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(3)), act_err_1(1,3), ['x$10^{' num2str(act_err_2(1,3)) '}$' ],  act_err_1(2,3), ['x$10^{' num2str(act_err_2(2,3)) '}$' ], act_err_1(3,3), ['x$10^{' num2str(act_err_2(3,3)) '}$' ], act_err_1(4,3), ['x$10^{' num2str(act_err_2(4,3)) '}$' ], act_err_1(5,3), ['x$10^{' num2str(act_err_2(5,3)) '}$' ], act_err_1(6,3), ['x$10^{' num2str(act_err_2(6,3)) '}$' ], act_err_1(7,3), ['x$10^{' num2str(act_err_2(7,3)) '}$' ]);
fprintf('$\\tau$ & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(4)), act_err_1(1,4), ['x$10^{' num2str(act_err_2(1,4)) '}$' ],  act_err_1(2,4), ['x$10^{' num2str(act_err_2(2,4)) '}$' ], act_err_1(3,4), ['x$10^{' num2str(act_err_2(3,4)) '}$' ], act_err_1(4,4), ['x$10^{' num2str(act_err_2(4,4)) '}$' ], act_err_1(5,4), ['x$10^{' num2str(act_err_2(5,4)) '}$' ], act_err_1(6,4), ['x$10^{' num2str(act_err_2(6,4)) '}$' ], act_err_1(7,4), ['x$10^{' num2str(act_err_2(7,4)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(5)), act_err_1(1,5), ['x$10^{' num2str(act_err_2(1,5)) '}$' ],  act_err_1(2,5), ['x$10^{' num2str(act_err_2(2,5)) '}$' ], act_err_1(3,5), ['x$10^{' num2str(act_err_2(3,5)) '}$' ], act_err_1(4,5), ['x$10^{' num2str(act_err_2(4,5)) '}$' ], act_err_1(5,5), ['x$10^{' num2str(act_err_2(5,5)) '}$' ], act_err_1(6,5), ['x$10^{' num2str(act_err_2(6,5)) '}$' ], act_err_1(7,5), ['x$10^{' num2str(act_err_2(7,5)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\\n', num2str(tau_vec(6)), act_err_1(1,6), ['x$10^{' num2str(act_err_2(1,6)) '}$' ],  act_err_1(2,6), ['x$10^{' num2str(act_err_2(2,6)) '}$' ], act_err_1(3,6), ['x$10^{' num2str(act_err_2(3,6)) '}$' ], act_err_1(4,6), ['x$10^{' num2str(act_err_2(4,6)) '}$' ], act_err_1(5,6), ['x$10^{' num2str(act_err_2(5,6)) '}$' ], act_err_1(6,6), ['x$10^{' num2str(act_err_2(6,6)) '}$' ], act_err_1(7,6), ['x$10^{' num2str(act_err_2(7,6)) '}$' ]);
fprintf('       & %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s & %1.2f %s \\\\ \\hline\n', num2str(tau_vec(7)), act_err_1(1,7), ['x$10^{' num2str(act_err_2(1,7)) '}$' ],  act_err_1(2,7), ['x$10^{' num2str(act_err_2(2,7)) '}$' ], act_err_1(3,7), ['x$10^{' num2str(act_err_2(3,7)) '}$' ], act_err_1(4,7), ['x$10^{' num2str(act_err_2(4,7)) '}$' ], act_err_1(5,7), ['x$10^{' num2str(act_err_2(5,7)) '}$' ], act_err_1(6,7), ['x$10^{' num2str(act_err_2(6,7)) '}$' ], act_err_1(7,7), ['x$10^{' num2str(act_err_2(7,7)) '}$' ]);
