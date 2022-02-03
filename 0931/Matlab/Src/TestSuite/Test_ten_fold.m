function Test_ten_fold()
% Compute the multiplicity structure for the ten fold system for n = 5 to 9 :
%     x_1 +...+ x_{i}  =  0, for i = 1, 2,..., n-2,
%     (x_1 + ...  + x_{n-1})^2  =  0, 
%     (x_1 + \cdots + x_{n})^5  =  0 at (0,...,0).
% The calling sequence is:
%     Test_ten_fold()

fid=fopen('ten_fold_time','wt');
disp('====================ten fold system===================');
for n=5:9
    disp(['====================','n = ',num2str(n),'===================']);
    [f,variables, zeros]=ten_fold(n);
    disp(['   system: ' f{1}])
    for j=2:length(f)
        disp(['           ' f{j}])
    end
    pause(1)
    disp(['   zero: ' num2str(zeros)])
    pause(1)
    disp(' ')
    disp('   Start computing the multiplicity structure')
    t1=cputime;
    [m,DI,HF]=multiplicity(f,variables,zeros);
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
