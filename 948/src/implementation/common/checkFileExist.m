function checkFileExist(filename)


fileStr = [pwd '/' filename];
err = exist(fileStr, 'file');

if err==2
    error('File ''%s'' exists', filename);
end

fileStr = [pwd '\' filename];
err = exist(fileStr, 'file');

if err==2
    error('File ''%s'' exists', filename);
end

end