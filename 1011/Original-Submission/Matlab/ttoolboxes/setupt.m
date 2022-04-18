function setupt( testflag, pathflag )
% setupt( [testflag], [pathflag] )
% Tests and install the packages: 
%       cell, double, m, sequence, subdivision, tjsr, tmisc, TTEST
%       The JSR Toolbox (JSR Louvain) v1.2b, SeDuMi 1.1,
% 
% Options:
%       testflag        0  ... fast test (default)
%                       1  ... slow test wi
%                       2  ... slow test without confirmation
%       pathflag        0  ... adds all packages to Matlab the path which do not seem to be installed (default)
%                       1  ... adds all packages to Matlab the path which do not seem to be installed and saves the matlab path
%                       2  ... asks whether to add packages to the Matlab path which do not seem to be installed and saves the Matlab path
%                       3  ... asks whether to add packages to the Matlab path and saves the Matlab path
%
% Written by: tommsch, 2019

%#ok<*NOPRT>
%#ok<*ASGLU>
%#ok<*NASGU>


    %parse input
    if( nargin<=0 || isempty(testflag) ); 
        testflag = 0; end;    
    if( nargin<=1 || isempty(pathflag) );
        pathflag = 0; end;
    
    fprintf( ['=====================================================================================\n', ...
              't-packages v1.0.7.1\n=====================================================================================\n'] );
          
	%pre-processing
    curdir = cd;
    if( exist(fullfile(pwd, 'setupt'), 'file')~=2 );
        cd( fileparts(which('setupt')) );
        fprintf( 'Warning: Function should be called from the directory where  ''setupt.m'' is located.\n' ); end;
    
    %set path
    addflag = false; 
    addflag = addflag | addpackage( 'EXPECT_EQ.m', '../TTEST', 'TTEST', pathflag ); 
    addflag = addflag | addpackage( 'cart2sphm2.m', 'm', 'm', pathflag ); 
    addflag = addflag | addpackage( 'symbol2mask.m', 'sequence', 'sequence', pathflag ); 
    addflag = addflag | addpackage( 'blf.m', 'subdivision', 'subdivision', pathflag ); 
    addflag = addflag | addpackage( 'tjsr.m', 'tjsr', 'tjsr', pathflag ); 
    addflag = addflag | addpackage( 'vprintf.m', 'tmisc', 'tmisc', pathflag ); 
    addflag = addflag | addpackage( 'sedumi.m', {'../sedumi','../sedumi/conversion','../sedumi/o_win'}, 'sedumi', pathflag );
    addflag = addflag | addpackage( 'jsr_pathcomplete.m', {'../JSR_louvain/Methods','../JSR_louvain','../JSR_louvain/Benchmark/','../JSR_louvain/Pre-processing','../JSR_louvain/Subroutines'}, 'JSR Louvain', pathflag );
    addflag = addflag | addpackage( pwd ); %special handling of the folder we are in, due to @cell, etc. packages
    
    
    if( ~testflag );
        fprintf( 'Matlab version: %s\n', version );
    else
        %get win/linux/mac
        fprintf( 'System architecture: %s\n', computer('arch') );
        
        %get processer
        if( ispc )      % Windows
            command = 'for /f "tokens=2 delims==" %A in (''wmic cpu get name /value'') do @(echo %A)';
        elseif( ismac ) % Mac
            command = 'sysctl machdep.cpu | grep brand_string | cut -d: -f2';
        else         % Linux
            command = 'grep -m 1 "model name" /proc/cpuinfo | cut -d: -f2'; end;
        [status, cpuName] = system( command );
        cpuName = strtrim( cpuName );
        if( status )
            cpuName = 'could not retrieve information'; end;
        fprintf( 'CPU name: %s\n', cpuName );
        
        % get matlab version and toolboxes
        ver; end;
    
    %post-pre-processing
    cd(curdir);
    
    
    %initialize TTESTs
    INIT( 'clear',1, 'all',testflag, 'v',0 );
    
    %print text
    fprintf( ['=====================================================================================\n'] );

    if( testflag>=1 )
        fprintf( ['This function runs a full test of all functions contained in the t-packages.\n', ...
                  'This test takes a long time. If it succeeds, you can be sure, that everything works as it should.\n'] );
        if( testflag<2 )
                  fprintf( 'Do you want to proceed (Press any key to continue, or ''Ctrl-C'' to abort) ?\n' ); end;
    else
        fprintf( ['This function runs a fast test of all functions contained in the t-packages.\n', ...
                  'To make a full test, call the function with parameter 1, i.e. ''setupt(1)''.\n'] ); end;

    if( isequal(testflag,1) );
        pause; end;
    fprintf( '=====================================================================================\n' );
    
    try; %#ok<TRYNC>
        evalc('plot([2]);'); %trigger possible warning: Warning: MATLAB has disabled some advanced graphics rendering features by switching to software OpenGL....
        close all; end;
    
    
    try
        %start tests    
        versionbad = test( 'testversion' );
        cellbad = test( 'testcell' );
        doublebad = test( 'testdouble' );
        mbad = test( 'testm' );
        sequencebad = test( 'testsequence' );
        subdivisionbad = test( 'testsubdivision' );
        tjsrbad = test( 'testtjsr' );
        tmiscbad = test( 'testtmisc' );

        bad = versionbad||cellbad||mbad||doublebad||sequencebad||subdivisionbad||tmiscbad||tjsrbad;
        fprintf( '=================\nSummary: ' );
        if( versionbad );     fprintf( '\n      Some toolboxes/external dependencies are not installed/licensed or Matlab version is not supported. Some parts of the t-toolboxes may not work.\n' ); end;
        if( cellbad );        fprintf( '\n      cell-package test run did not succeed.\n' );                  end;
        if( doublebad );      fprintf( '\n      double-package test run did not succeed.\n' );                end;
        if( mbad );           fprintf( '\n      m-package test run did not succeed.\n' );                     end;
        if( sequencebad );    fprintf( '\n      sequence-package test run did not succeed.\n' );              end;
        if( subdivisionbad ); fprintf( '\n      subdivision-package test run did not succeed.\n' );           end;
        if( tjsrbad );        fprintf( '\n      tjsr-package test run did not succeed.\n' );                  end;
        if( tmiscbad );       fprintf( '\n      tmisc-package test run did not succeed.\n' );                 end;
        if( ~bad );           fprintf( 'If there are no errors/warning the test succeeded.\n' );              end;
    catch
        fprintf( 'Fatal error while testing.' ); end;
    
    if( pathflag>=1 );
        savepath; 
    elseif( addflag )
        fprintf( '\nTo use the t-packages regulary, \nyou must save this new path definition.\n');
        fprintf( 'To do this, type the command \n     ''savepath'' \nat the Matlab prompt. \nPlease consult the MATLAB \ndocumentation if necessary.\n'); end;
    
   
    close all;
end

function flag = addpackage( file, relpath, name, pathflag )
% flag is set, if something is added to matlab path

flag = false;

    if( nargin==4)
    

        if( ~iscell(relpath) );
            relpath = {relpath}; end;

        if( exist(fullfile(pwd, 'setupt'), 'file')~=2 );
            error( 'Function must be called from the directory where  ''setupt.m'' is located.' ); end;

        if( exist(file,'file')~=2 || pathflag>=3 )
            if( pathflag>=2 );
                x = '';
                while( ~isequal(x,'y') && ~isequal(x,'n') );
                    x = input( ['Do you want to add the'  name 'package to the path (y/n)?'], 's' ); end;
            else;
                x = 'y'; end;
            if( isequal(x,'y') );
                if( isequal(exist(fullfile(relpath{1}, file), 'file'),2) );
                    for i = 1:numel( relpath );
                        p = fullfile( pwd, relpath{i} );
                        fprintf( 'Add to matlab-path: %s\n', p );
                        flag = true;
                        addpath( p ); end;
                else
                    fprintf( 'ERROR: Missing files. Package %s not added to path.\n', name ); end; end; 
        else
            %fprintf( 'Package already installed: %s\n', name ); 
            end;
    elseif( nargin==1 )
        if ispc  % Windows is not case-sensitive
          %onpath = any(strcmpi(file, str));
          onpath = ~isempty( strfind(lower(path),[lower(file) pathsep]) ); %#ok<STREMP>
        else
          %onpath = any(strcmp(file, pathCell)); 
          onpath = ~isempty( strfind(path,[file pathsep]) ); end; %#ok<STREMP>
        if( ~onpath );
            flag = true;
            addpath( file ); end; end;

end

function ret = test( str )
    ret = runtests(str);
    ret = ~all([ret.Passed]);
    if( flag )
        if( INIT('f') );
            ret = true;
            fprintf( 'Test failed\n' ); end; end;
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 
