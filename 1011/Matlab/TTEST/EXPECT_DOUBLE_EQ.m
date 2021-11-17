%TTEST_AUTOGENERATE
%This file is auto-generated. Do not modify it. Changes may be overwritten.
function ret = EXPECT_DOUBLE_EQ( varargin )

global TTEST_ERROR_EXPECT;
global TTEST_ERRORFLAG;

ret = true;
try;
    if( nargin==2 );
        factor = 4; 
    else
        factor = varargin{3}; end;

    val = norm( varargin{1}-varargin{2} );
    maxdiff = factor*eps(double(norm(varargin{1})));
    if( val>maxdiff );
        ret = false;
        TTEST_FPRINTF( '\n\nlhs != rhs (but they should be almost equal), where');
        TTEST_FPRINTF( '\nlhs = ' );
        TTEST_DISP( varargin{1} );
        TTEST_FPRINTF( '\nrhs = ' );
        TTEST_DISP( varargin{2} );
        TTEST_FPRINTF( '\nDifference norm(lhs-rhs) = %.15g\n', val );
        TTEST_FPRINTF( 'Allowed difference       = %.15g', maxdiff ); end;
catch;
    error( 'TTEST_EVAL:format', 'Wrong format of varargin' ); end;

if( isequal(ret,false) );
    TTEST_ERRORFLAG = true;

    try;
        if( isequal(TTEST_ERROR_EXPECT,'exp') );
            TTEST_ERROR_EXPECT_WORKER();
        elseif( isequal(TTEST_ERROR_EXPECT,'ass') );
            TTEST_ERROR_ASSERT_WORKER();
        else;
            feval( TTEST_ERROR_EXPECT ); end;
    catch me;
        if( isequal(me.identifier,'MATLAB:assertion:failed') );
            assert( false ); 
        else; 
            rethrow( me ); end; end; end;
end

