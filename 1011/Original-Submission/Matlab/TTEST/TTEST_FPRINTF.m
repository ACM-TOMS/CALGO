function TTEST_FPRINTF( varargin );
    global TTEST_VERBOSE;
    if( isempty(TTEST_VERBOSE) || TTEST_VERBOSE>= 1 );
        fprintf(varargin{:}); end;
end