% append_MinGW_dir.m
%{
% Under Windows:
%     this script appends, only temporarily, the MinGW binary folder (MinGW/bin)
%     to the system PATH variable; "temporarily" means that, at MATLAB exit,
%     the system PATH variable returns as before.
%     This is required to run C/MATLAB mixed code (otherwise an "Invalid
%     mex file" error is issued).
%     If C uses the OpenMP directives then this folder must be added
%     permanently.
%}

fprintf('\n****************************************************************************************\n')
fprintf('                               append_MinGW_dir utility\n')
if ispc % ONLY if is Windows
    fprintf('\nunder Windows:\n')
    PATHnow = getenv('PATH');
    if isempty(strfind(PATHnow, '\mingw64\bin')) % if MinGW gcc compiler's folder not yet added
        fprintf('the directory \\mingw64\\bin is not in the PATH environment variable.\n')
        fprintf('Now search in C:\\Program Files\n'); % default directory
        lst=dir('C:\Program Files');
        notFound=isempty(strfind(lst(1).name,'mingw-w64'));
        if notFound % search
            k=2;
            while k<=numel(lst)
                notFound=isempty(strfind(lst(k).name,'mingw-w64'));
                if ~notFound % found!
                    break
                end
                k=k+1;
            end
        end
        if k > numel(lst)
            disp('NOT FOUND!')
            answer = inputdlg('Enter the complete path of bin subfolder for MinGW installation:','MinGW',[1 80], ...
                              {'C:\Program Files\mingw-w64\x86_64-5.1.0-posix-seh-rt_v4-rev0\mingw64\bin'});
            if isempty(answer)
                error('*** Entered an empty folder. Install before MinGW! ***')
            else
                disp('Now execute:')
                disp(['>> setenv(PATHnow;' answer '\mingw64\bin)'])
                setenv('PATH', [PATHnow ';C:\Program Files\mingw-w64\' lst(k).name '\mingw64\bin']);
                fprintf('Appended to the PATH environment variable the directory /bin of MinGW 64bit\n')
            end
        else
            disp(['"mingw-w64" FOUND in lst(' num2str(k) ')'])
            lst=dir('C:\Program Files\mingw-w64');
            for k=1:numel(lst)
                if ~strcmp(lst(k).name,'.')  &&  ~strcmp(lst(k).name,'..')
                    break
                end
            end
            disp('Now execute:')
            disp(['>> setenv(PATHnow;C:\Program Files\mingw-w64\' lst(k).name '\mingw64\bin)'])
            setenv('PATH', [PATHnow ';C:\Program Files\mingw-w64\' lst(k).name '\mingw64\bin']);
            fprintf('Appended to the PATH environment variable the directory /bin of MinGW 64bit\n')
        end
    else
        fprintf('the directory \\mingw64\\bin already in the PATH environment variable.\n')
    end
else
    fprintf('\nMATLAB is not running under Windows.\n')
end
fprintf('\n****************************************************************************************\n')
