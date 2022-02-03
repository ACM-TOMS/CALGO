% DISP    DISP(OBJ) displays object obj in MATLAB syntax
function disp(obj)
    c = char(obj);
    if iscell(c)
        disp(['     ' c{:}])
    else
        disp(c)
    end
end
