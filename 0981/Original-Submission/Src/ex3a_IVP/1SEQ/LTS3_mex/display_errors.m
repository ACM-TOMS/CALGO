% display_errors.m

%% DISPLAY ERROR MATRICES
disp(['RELERR1 = [ %  Tval(1)    Tval(2)    ...    Tval(' num2str(NTval) ')']);
disp([repmat(' ',NXval,11) num2str(RELERR1,'%13.6e') repmat('    % Xval(',NXval,1) num2str((1:NXval)','%5d') repmat(')',NXval,1)]); fprintf('          ];\n')
fprintf('\n')
disp(['RELERR2 = [ %  Tval(1)    Tval(2)    ...    Tval(' num2str(NTval) ')']);
disp([repmat(' ',NXval,11) num2str(RELERR2,'%13.6e') repmat('    % Xval(',NXval,1) num2str((1:NXval)','%5d') repmat(')',NXval,1)]); fprintf('          ];\n')
