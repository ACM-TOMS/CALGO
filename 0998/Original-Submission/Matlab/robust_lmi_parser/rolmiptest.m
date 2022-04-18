clear;

try
    P = sdpvar(2,2,'sym');
    clear P;
    P = rolmipvar(2,2,'P','sym',2,1);
    clear P;
    
    disp('ROLMIP is properly installed.');
catch ME
    if (strcmp(ME.message,'Undefined function ''sdpvar'' for input arguments of type ''double''.'))
        fprintf(2,'YALMIP Toolbox not installed! In order to properly run ROLMIP, \n it is necessary to install YALMIP, which can be downloaded at \n https://yalmip.github.io/download/ \n \n');
    elseif (strcmp(ME.message,'Undefined function ''rolmipvar'' for input arguments of type ''double''.'))
        fprintf(2,'ROLMIP is not properly installed! Please check if the \n installation directory is included into MATLAB path. \n \n');
    end
end
        


