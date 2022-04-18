function Test_cyclic_cubic()
% Compute the multiplicity structure for the cyclic cubic system for n = 5 to 7 :
%    x^3_i - x_{i+1} x_{i+2} = 0, for i = 1, 2,..., n, 
%    x_{n+1} = x_1 and x_{n+2} = x_2 at zero (0,...,0).
% The calling sequence is:
%     Test_cyclic_cubic()

fid=fopen('cyclic_cubic_time','wt');
disp('====================Cyclic cubic system===================');
for n=5:7
    disp(['====================','n = ',num2str(n),'===================']);
    [f,variables,zeros,options]=cyclic_cubic(n);
    disp(['   system: ' f{1}])
    for j=2:length(f)
        disp(['           ' f{j}])
    end
    pause(1)
    disp(['   zero: ' num2str(zeros)])
    pause(1)
    disp(' ');
    disp('   Start computing the multiplicity structure')
    t1=cputime;
    [m,DI,HF]=multiplicity(f,variables,zeros,options);
    t2=cputime;
    disp(['   computation finished in ',num2str(t2-t1),' seconds, with following results:']);
    disp(' ');
    pause(1)
    disp(['   multiplicity: ' num2str(m)])
    pause(1)
    disp('   the dual basis: ')
    outputDI(DI)
    pause(1)
    disp(['   the Hilbert function: ' num2str(HF)])    
    disp(' ')
end
fclose(fid);
