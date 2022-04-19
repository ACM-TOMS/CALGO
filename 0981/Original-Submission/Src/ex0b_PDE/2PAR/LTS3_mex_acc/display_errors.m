% display_errors.m

%% DISPLAY ERROR MATRICES
if ABSERR
    ERRstr='ABSERR';
else
    ERRstr='RELERR';
end

fprintf('\n')
disp([ERRstr '1 = [ %  Tval(1)     Tval(2)      Tval(3)      Tval(4)  ... Tval(' num2str(NTval) ')']);
disp([repmat(' ',NXval,11) num2str(ERROR1,'%13.6e') repmat('    % Xval(',NXval,1) num2str((1:NXval)','%5d') repmat(')',NXval,1)]); fprintf('          ];\n')

fprintf('\n')
disp([ERRstr '2 = [ %  Tval(1)     Tval(2)      Tval(3)      Tval(4)  ... Tval(' num2str(NTval) ')']);
disp([repmat(' ',NXval,11) num2str(ERROR2,'%13.6e') repmat('    % Xval(',NXval,1) num2str((1:NXval)','%5d') repmat(')',NXval,1)]); fprintf('          ];\n')

fprintf('\n')
disp([ERRstr '3 = [ %  Tval(1)     Tval(2)      Tval(3)      Tval(4)  ... Tval(' num2str(NTval) ')']);
disp([repmat(' ',NXval,11) num2str(ERROR3,'%13.6e') repmat('    % Xval(',NXval,1) num2str((1:NXval)','%5d') repmat(')',NXval,1)]); fprintf('          ];\n')
