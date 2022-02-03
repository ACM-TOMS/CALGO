% display_time_perc.m
%% time percentages: time/tot_time

%% DISPLAY TIME PERCENTAGE MATRICES
disp('PARperc1 = [ %            5           20          120  = NXval');
disp([repmat(' ',numel(NTval),15) num2str(PARtime1./TOTtime1,'%13.6e') repmat('  % ',numel(NTval),1) num2str(NTval,'%5d') repmat(' = NTval',numel(NTval),1)]); fprintf('            ];\n')
disp('LTSperc1 = [ %            5           20          120  = NXval');
disp([repmat(' ',numel(NTval),15) num2str(LTStime1./TOTtime1,'%13.6e') repmat('  % ',numel(NTval),1) num2str(NTval,'%5d') repmat(' = NTval',numel(NTval),1)]); fprintf('            ];\n')
disp('SUMperc1 = [ %            5           20          120  = NXval');
disp([repmat(' ',numel(NTval),15) num2str(SUMtime1./TOTtime1,'%13.6e') repmat('  % ',numel(NTval),1) num2str(NTval,'%5d') repmat(' = NTval',numel(NTval),1)]); fprintf('            ];\n')

disp('PARperc2 = [ %            5           20          120  = NXval');
disp([repmat(' ',numel(NTval),15) num2str(PARtime2./TOTtime2,'%13.6e') repmat('  % ',numel(NTval),1) num2str(NTval,'%5d') repmat(' = NTval',numel(NTval),1)]); fprintf('            ];\n')
disp('LTSperc2 = [ %            5           20          120  = NXval');
disp([repmat(' ',numel(NTval),15) num2str(LTStime2./TOTtime2,'%13.6e') repmat('  % ',numel(NTval),1) num2str(NTval,'%5d') repmat(' = NTval',numel(NTval),1)]); fprintf('            ];\n')
disp('SUMperc2 = [ %            5           20          120  = NXval');
disp([repmat(' ',numel(NTval),15) num2str(SUMtime2./TOTtime2,'%13.6e') repmat('  % ',numel(NTval),1) num2str(NTval,'%5d') repmat(' = NTval',numel(NTval),1)]); fprintf('            ];\n')
