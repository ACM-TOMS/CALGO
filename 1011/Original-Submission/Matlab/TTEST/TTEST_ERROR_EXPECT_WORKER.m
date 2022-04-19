function TTEST_ERROR_EXPECT_WORKER
% TTEST_ERROR_EXPECT_WORKER()
% Default function called when an EXPECT_ fails
%
% Tries to print out the filename, linenumber and testname
%
% See also: TTEST_ERROR_ASSERT_WORKER

    ds = dbstack;
    stackidx = min(3, numel(ds) );
    
    name =  ds(stackidx).file;
    fullname = which(name);
    fprintf( '\nIn: <a href="matlab: opentoline(''%s'',%i,0)"  style="font-weight:bold">%s:%i</a>\n', fullname, ds(stackidx).line, name, ds(stackidx).line );
    str = get_current_section_name;
    if( ~isempty(str) );
        fprintf( 'Testname: %s\n', get_current_section_name ); end;
    fprintf( '=====================================\n' );
end


function sectionname = get_current_section_name

%get string of section name

ds = dbstack();
stackidx = min(4, numel(ds) );

execution_file = ds(stackidx).file;
execution_line = ds(stackidx).line;

fid = fopen( execution_file );
sectionname = '';
current_line = 1;

while( true)
    tline = fgetl( fid );
    if( ~ischar(tline) );
        sectionname = '';
        break; end;
    if( startsWith(tline,'%%') )
        sectionname = tline; end;
    if( execution_line <= current_line )
        break; end;
    current_line = current_line + 1; end;

fclose( fid );

if( numel(sectionname)>=3 );
    sectionname = sectionname(3:end); end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 