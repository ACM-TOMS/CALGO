function Test_KSS()
% Compute the multiplicity structure for the KSS system for n = 5 to 8 :
%     x_i^2+\sum_{j=1}^{n}x_j-2x_i-(n-1) = 0, for  i=1,2,...,n at the zero (1,...,1).
% The calling sequence is:
%     Test_KSS()

fid=fopen('KSS_time','wt');
disp('====================KSS system===================');
for n=5:8
    disp(['====================','n = ',num2str(n),'===================']);
    [f,variables,zeros,options]=KSS(n);
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
