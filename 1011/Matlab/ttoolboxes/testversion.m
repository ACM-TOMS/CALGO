% test version

%#ok<*CTPCT>
%#ok<*ASGLU>
%#ok<*NASGU>

% preconditions    

%% test external functions / libraries

%% test version
assert( ~verLessThan('matlab','9.1') ); 

%% test parallel toolbox

inst = ~isempty( ver('distcomp') );
lic = license( 'test', 'Distrib_Computing_Toolbox' );
file = isequal( exist('parfor','file'), 2 );
if( ~lic || ~file || ~inst );
    fprintf('\n');
    warning( 'testtjsr:distcomp', 'ERROR: Parallel Toolbox is not installed or licensed. Some functions in the t-toolbox may not work.' ); 
    assert(false);
else
    gcp(); %start parallel pool
end;

    
%% test symbolic toolbox
inst = ~isempty( ver('symbolic') );
lic = license('test', 'symbolic_toolbox');
file = isequal( exist('sym','file'), 2 );
if( ~lic || ~file || ~inst );
    fprintf('\n');
    warning( 'testsequence:symbolic', 'ERROR: Symbolic Toolbox is not installed or licensed. Some functions in the t-toolbox may not work.' ); 
    assert(false); end;

%% test signal processing toolbox
inst = ~isempty( ver('signal') );
lic = license('test', 'signal_toolbox');
file = isequal( exist('parfor','file'), 2 );
if( ~lic || ~file || ~inst );
    fprintf('\n');
    warning( 'testsequence:signal', 'ERROR: Signal Processing Toolbox is not installed or licensed. Some functions in the t-toolbox may not work.' ); 
    assert(false); end;



%% test gurobi
try
    evalc( 'gurobi_setup' ); % run evalc to supress output
    model.A = sparse( [1 1 0; 0 1 1] );
    model.obj = [1 2 3];
    model.modelsense = 'Max';
    model.rhs = [1 1];
    model.sense = [ '<' '<' ];
    params.outputflag = 0;

    gurobiresult = gurobi( model, params );
    if( ~isequal(gurobiresult.x,[ 1 0 1].') ); 
        fprintf( 'ERROR: Gurobi broken, installed but not working.\n' );
        fprintf( 'ERROR: ==========================\nThe tjsr-package will NOT work and may NOT emit error-messages.\n' );
        assert( false ); end; 
catch
    fprintf('\n');
    fprintf( 'Gurobi is not installed. tjsr will use the Matlab solver linprog.\n' );
    %warning( 'testversion:gurobi', 'Gurobi is not installed. The tjsr toolbox will work, but magnitudes slower.' );
    %assert(false); 
    end;
   
%% test JSR louvain
try
    A = {[1 -1 0;2 -1 -1; 0 1 -1],[1 3 -4;2 -3 1;4 -2 -2]};
    evalc( 'b = jsr_pathcomplete(A);' );
    b = (b > 3.605 && b < 3.61);
    a = jointTriangul( A );
    
    assert( a&&b );
catch
    warning( 'testversion:jsrlouvain', 'The JSR Toolbox (JSR Louvain) is not installed/working. Some functions in the tjsr toolbox may not work.\n' );
    assert( false ); end;


%% test sedumi
try
    AA = [ 4  4  4 -1 -4 -4
           4  1 -4 -4 -4  4]';
    bb = [ 2  1]';
    cc = [-4  0  0  0  0  0]';
    pars.fid = 0;
    K.l = 6;
    val = sedumi( AA, bb, cc, K, pars );
    assert( norm( val - [0.250000000000000   0.200025626486528   0.178778846360113   0.211428602830903   0.242873139784766   0.333074182354149]' )<1e-12 ); 
catch;
    warning( 'testversion:sedumi', 'SeDuMi not installed. Some functions in the tjsr toolbox may not work.\n' );
    assert( false ); end;
