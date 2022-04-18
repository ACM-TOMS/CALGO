function msg(logfile,printCond,string,varargin)
%
% MSG(LOGFILE,PRINTCOND,STRING,VARARGIN)
%    will show the message STRING with eventual parameters in VARARGING
%    on the command line if PRINTCOND (boolean). And will always write it
%    in the file corresponding to the file identifier LOGFILE created by
%    fopen.
%      
    
if (printCond)
disp(sprintf(string, varargin{:}));
end  

if (logfile ~= -1)
    fprintf(logfile,sprintf([ string '\n'], varargin{:}));
end