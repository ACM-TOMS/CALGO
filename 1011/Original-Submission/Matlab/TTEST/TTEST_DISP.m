function TTEST_DISP( x )
    global TTEST_VERBOSE;
    global TTEST_MAXDISP;
    if( isempty(TTEST_VERBOSE) || TTEST_VERBOSE >= 1);
        str = strtrim(evalc('disp(x)'));
        if( isempty(x) );
            if( iscell(x) );
                str = '{}';
            else;
                str = '[]'; end;
        elseif( isempty(TTEST_MAXDISP) || numel(str) > TTEST_MAXDISP );
            str(TTEST_MAXDISP+1:end) = [];
            str = [str ' ... Output truncated' newline]; end;
    fprintf('%s',str); end;
end