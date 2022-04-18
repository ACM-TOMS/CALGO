%% Create Source Code Documentation
% Process to create HTML files for source code documentation in the 
% directory |doc/html| out of in-place documentation in MATLAB files of the 
% toolbox |IPscatt| by the function |publish| of MATLAB. This routine can
% be used for other toolboxes too.
%
% Warning: This function creates a directory and files.
%
%% Syntax
%
%   docCreate;
%
%% Description
%
% The routine |docCreate| uses the function |publish| of MATLAB to publish
% the source code documentation written as in-place documentation of the
% toolbox |IPscatt|.
% 
% In detail:
%
% * Takes into account all m-files in the folder |code| and its subfolders.
% * Uses the function |publish| of MATLAB to create HTML-files in the
% folder |doc/html|.
% * Creates an index of all routines, |index.html|, in the folder |doc|.
% * Takes into account additional pictures stored in |doc/extGraph| when referred in the documentation.
% Actually, these is an overview of the program's structure and a 
% reference of all files (created manually copying all input and output
% arguments of the MATLAB files of |IPscatt|.)
% * Includes additional MATLAB-files in |code/docCreate/addMfiles|.
% (The files |programStructure| and |setiRef| have to be updated manually,
% if the structure in form of file dependencies is changed or references
% of input and output arguments are changed.)
%
% To save time in case of updates of the source code documentation the
% files in the folder |html| and the file |index.html| are not deleted,
% such that it may contain documentation of deleted files. 
% They can be deleted without data loss.
%
%% Example
%
% *Example 1: Generate or update the documentation of all files and create the index file*
%
%   docCreate;
%
% *Example 2: Generate or update the documentation of specific file*
%
% This example updates the source code documentation of the file |start.m|.
%
%   options.outputDir = '../doc/html';
%   options.format    = 'html'; % default
%   options.evalCode  = false;  % Do not evaluate the code.
%   publish('start.m',options);
%
%% More About
%
% *Structure*: In general we follow the structure of the MATLAB documentation. i.e.
%
% # Name of the function and a short description.
% # Syntax
% # Description
% # Examples
% # Input Arguments
% # Output Arguments
% # More About
% # References
% # See Also
%
% Additionaly we have the section Code.
%
% *URLs*: Some interesting URLs how to write an in-place documentation in MATLAB:
%
% * <http://de.mathworks.com/help/matlab/matlab_prog/publishing-matlab-code.html>
% * <http://de.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html>
% * <http://de.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html>
% * <https://wiki.hpcc.msu.edu/pages/viewpage.action?pageId=17599292>
% (auto-documentation in MATLAB)
%
% *Others*:
%
% * Currently, all folders with the string |closed| are ignored for the 
% publication of the source code documentation. A distinction between |doc| 
% and |docClosed| in dependence of |closed = 0| or |closed = 1| would be
% useful, but is not implemented yet. This should take into account the
% added paths in paths in |init.m|.
%
% *Look for errors*
%
% LaTeX may fail, e.g.
%
%   $$\alpha = \alpa$$
%
% Instead of checking each file, look for "Error updating Text." in HTML
% files, e.g. with "Find Files" from MATLAB.
%
%% Code: docCreate

function docCreate()

init;

if exist('closed','var') && closed == 1
    disp('Parameter closed is 1 in init.m.')
    disp('The source code documentation was only written for the public version.')
    error('Do not set the variable closed in init or set them to 0 and run docCreate again.')
end

pubMfiles = 1; % start publishMfiles (takes some time...)

%% 
% *What happens?*
%
% * Make make a link for each file by
%
%   <a href='...'>...</a>
%
% * Structure them (folder-wise).
%

%%
% *Create a list of all folders*

disp(' ')
disp('## List of all m-files --------------------------------------------')
disp(' ')

%%% Look for folders
disp(' ')
disp('-- List of folders --')
disp(' ')
[folderlist,~] = createFolderList();
% transpose(folderlist)

disp(' ')
disp('-- Skip some items in list of folders --')
disp(' ')
folderlist = skipItemsInFolderList(folderlist);
% transpose(folderlist)

%%% Look for m-files
disp(' ')
disp('-- List of m-files --')
disp(' ')

[mfileDirname,mfileName,mfilePath] = createMfileList(folderlist);

%%
% *Make index*

disp(' ')
disp('## Make index HTML file -------------------------------------------')
disp(' ')
makeIndexHTML(mfileDirname,mfileName);

%%
% *Publish it*

disp(' ')
disp('## Publish it -----------------------------------------------------')
disp(' ')
if pubMfiles
    options = pubOptions();
    publishMfiles(mfilePath,options);
end

end

%% Code: Subfunction: publishMfiles

function publishMfiles(mfilePath,options)
for i = 1:length(mfilePath)
    fprintf('   Publish %s',mfilePath{i})
    publish(mfilePath{i},options);
    fprintf(' -- [OK]\n')
end
end

%% Code: Subfunction: makeIndexHTML

function makeIndexHTML(mfileDirname,mfileName)

fid = fopen('../doc/index.html', 'wt');
%  fprintf( fid, '%f,%f,%f,%f\n', a1, a2, a3, a4);
fprintf(fid,'<html>\n');
fprintf(fid,'<head>\n');
fprintf(fid,'<title>Documentation: IPscatt</title>\n');

%%% Copy CSS from MATLAB Standard Docu into index.html
fprintf(fid,'<!-- css formats are copied from MATLAB published documentation -->\n');
fprintf(fid,'<style type="text/css">\n');
fidcss = fopen('docCreate/docCreateCSS.css','r');
tline = fgetl(fidcss);
while ischar(tline)
    %disp(tline)
    tline = fgetl(fidcss);
    fprintf(fid,'%s\n',tline);
end
fclose(fidcss);
fprintf(fid,'</style>\n');
%fprintf(fid,css);
%clear css;

fprintf(fid,'</head>\n');
fprintf(fid,'<body>\n');
fprintf(fid,'<div class="content">\n');
fprintf(fid,'<h1>Source Code Documentation of IPscatt&mdash;a MATLAB Toolbox for the Inverse Medium Problem in Scattering</h1>\n');
% fprintf(fid,'<br>\n');

fprintf(fid,'<p>This is the source code documentation of the toolbox IPscatt.</p>\n');

fprintf(fid,'<p>It is structured as follows:</p>\n');

fprintf(fid,'<ul>');

fprintf(fid,'<li>Read ');
htmlfile = 'start.html';
fprintf(fid,'<a href="html/%s">%s</a>',htmlfile,htmlfile);
fprintf(fid,' for a brief introduction in the usage of IPscatt and the underlying theory of direct and inverse scattering from inhomogeneous media.</li>\n');

fprintf(fid,'<li>See ');
htmlfile = 'programStructure.html';
fprintf(fid,'<a href="html/%s">%s</a>',htmlfile,htmlfile);
fprintf(fid,' for dependencies of functions.</li>\n');

fprintf(fid,'<li>Below is the structure of files in folders. \n');
fprintf(fid,'We provide a source code documentation for each function.</li>\n');

fprintf(fid,'<li>Furthermore, in ');
htmlfile = 'setiRef.html';
fprintf(fid,'<a href="html/%s">%s</a>',htmlfile,htmlfile);
fprintf(fid,' are fields of the main structure array seti summarized in one file. ');
fprintf(fid,'(The fields are used as input and output and described in the source code documentation of each corresponding file too.)');
fprintf(fid,'</li>\n');
fprintf(fid,'</ul>\n');

fprintf(fid,'<h2>Structure of Files in Folders of IPscatt</h2>');

dirnamePrev = 0;
for i = 1:length(mfileName)

    % Display foldername if new folder
    dirname = mfileDirname{i};
    if ~strcmp(dirnamePrev,dirname)
        fprintf(fid,'<br><h5>%s</h5>',dirname);
    end
    dirnamePrev = dirname;
    
    % Link to documentation file
    htmlFilename = strrep(mfileName{i},'.m','.html'); % Replace .m by .html for the path.
    fprintf(fid,'<a href="html/%s">%s</a><br>',htmlFilename,mfileName{i});
end

fprintf(fid,'</div>\n');
fprintf(fid,'</body>\n');
fprintf(fid,'</html>\n');
fclose(fid);

end

%% Code: Subfunction: createMfileList

function [listDirname,listName,listPath] = createMfileList(folderlist)
%%
% Output:
%
% * |listName|: name of m-file
% * |listPath|: path of m-file

display = 0; % 1: display results, otherwise 0
k = 0;

for i = 1:length(folderlist)+1
    if display == 1
        disp(' ')
    end
    if i == 1
        % files in current folder
        dirname = '';
        dirstring = '*.m'; % current folder
        mfile = dir(dirstring);
        if display == 1
            disp('*.m');
        end
        lookInFolder = 1;
    elseif ~isempty(folderlist{i-1}) % not used folders are marked by emty cell arrays
        % folderlist is from 1:length(folderlist)
        % using i-1 because i = 1 is reserved for first level
        dirname = folderlist{i-1};
        dirstring = sprintf('%s*.m',folderlist{i-1});
        mfile = dir(dirstring);
        if display == 1
            disp(dirstring);
        end
        lookInFolder = 1;
    else
        lookInFolder = 0;
    end
    
    if lookInFolder == 1
        for j = 1:length(mfile)
            mfileentry = sprintf('%s',mfile(j).name);
            %if exist(mfileentry,'file') == 2
            % Do not use the line above because
            % * procpost/planeplot.m and
            % * procpost/plot3Drewrite.m
            % are excluded (do not know why...)
                k = k + 1;
                if display == 1
                    fprintf('%03d | %02d | %s\n',k,j,mfileentry);
                end
                listDirname{k} = dirname; % store dirname
                listName{k} = mfileentry;
                listPath{k} = sprintf('%s%s',dirname,mfileentry);
            %end
            clear mfileentry;
        end
    end
end

end

%% Code: Subfunction: createFolderList

function [list,j] = createFolderList(varargin)
%%
% *Create a folder list*
%
% Store the relative path of (approximately) all folders and subfolders
%
%   [folderlist,j] = createFolderList();
% or
%   [folderlist,j] = createFolderList(list,j,path);
%   [folderlist,j] = createFolderList('',0,'proc/');

%%
% But this does not contain all folders..., e.g. not
%
% * inexpdata/fresnel_opus_1,
% * inexpdata/fresnel_opus_2.
%
% It is not clear why, but these folders are uninteresting for the
% documentation.

display = 0; % 1: display results, 0: not...

switch nargin
    case 0
        list = '';
        j = 0;
        folder = dir;
        path = '';
    case 1 || 2
        error('folderlist needs 0 or 3 arguments.')
    case 3
        list = varargin{1};
        j    = varargin{2};
        path = varargin{3};
        folder = dir(path); %e.g. folder = dir('proc')
end

for i = 1:length(folder)
    folderentry = sprintf('%s',folder(i).name);
    
    % Do not use
    % exist(folderentry,'file') == 7
    % instead of
    % exist(folderentry,'dir') == 7
    % Why?
    % 'file' does not work correctly...
    % Some dirs are not recognized, if you put the name as string.
    % (e.g. proc/setData).
    % (It is no problem in case of exist('proc/setData','file') == 7.)
    
    if exist(folderentry,'dir') == 7 && ~strcmp(folderentry,'.') && ~strcmp(folderentry,'..')
        j = j + 1;
        %list{j} = folderentry;
        if display == 1
            fprintf('%02d | %s%s\n',j,path,folderentry)
        end
        list{j} = sprintf('%s%s/',path,folderentry);
        [list,j] = createFolderList(list,j,list{j});
%     else
%         fprintf('NO FOLDER: %s%s\n',path,folderentry)
    end
    
    clear folderentry;
end
end

%% Code: Subfunction: skipItemsInFolderList

function folderlist = skipItemsInFolderList(folderlist)
for i = 1:length(folderlist)
    % Ignore folders with string "closed" by setting the item to [].
    str = folderlist{i};
    expression = 'closed'; % \w = [a-zA-Z_0-9]
    matchStr = regexp(str,expression,'match','ignorecase');
    if ~isempty(matchStr)
        folderlist{i} = []; % delete cell array entry
    end
    clear str expression matchStr;
    % ignore folder nogit
    if strcmp(folderlist{i},'nogit/')
        folderlist{i} = [];
    end
end
end

%% Code: Subfunction: pubOptions

function options = pubOptions()
p = mfilename('fullpath'); % full path of this file
f = fileparts(p); % folder path (without file)
pubDir = sprintf('%s/../doc/html',f); % full path
clear p f;

if exist(pubDir,'dir') ~= 7
    mkdir(pubDir);
end

%options.outputDir = '/home/fbuergel/Documents/rebis/doc/html'; % full path
options.outputDir = pubDir; clear pubDir;
options.format    = 'html'; % default
options.evalCode  = false; % do not evaluate the code

%options.showCode = true; % show code (default)
%options.showCode = false; % do not show code

%publish('start.m',options)
%clear options
end

%% Remember the function exist
% Remember the MATLAB function exist, see
% <https://de.mathworks.com/help/matlab/ref/exist.html>.
% 
% *2, i.e. one of the following is true*
%
% * name exists on your MATLAB search path as a file with extension .m or .mlx.
% * name is the name of an ordinary file on your MATLAB search path.
% * name is the full pathname to any file.
%
% *7, i.e. the name is a folder*
%
%   A = exist('name','kind')
%
% * dir : Checks only for folders.
% * file: Checks only for files or folders.

% if exist('...','file')
%     if exist('...','dir') % must be a folder
%        ...
%            
%     else % must be a file
%         dir *my*.m
%         
%     end
% end
% 
% 
