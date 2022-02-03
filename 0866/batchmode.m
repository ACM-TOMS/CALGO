function batchmode(testproblem)
% batchmode    enables batch processing for IFISS testproblem 
%   batchmode(testproblem);
%   input
%          testproblem  character string naming the testproblem
%                       must have the form "*_batch".m where "*" begins with
%                       "P"       for Poisson problems
%                       "CD"      for convection-diffusion problems
%                       "S"       for Stokes problems
%                       "NS"      for Navier-Stokes problems
%                       "itsolve" for iterative solution of linear systems
%   side effect
%          If batchmode terminates prematurely because of an error or
%          execution of cntl-C, interactive input with IFISS may not
%          work correctly.  This is fixed by typing "activemode".
%
%
%   IFISS function: DJS, HCE; 28 September 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 

global BATCH FID

% file containing input data 
batchfile=[testproblem,'_batch.m'];
[FID,message]=fopen(batchfile,'r');

% Error return on nonexistent or misnamed input file
if strcmp(message,'')~=1
   error(['INPUT FILE ERROR: ' message])
else
   disp(['Working in batch mode from data file ' batchfile])
end
if ~( strncmp(testproblem,'itsolve',7) | ...
      strncmp(testproblem,'P',1) | strncmp(testproblem,'CD',2) | ...
      strncmp(testproblem,'S',1) | strncmp(testproblem,'NS',2) ),
    errmsg = 'INPUT FILE ERROR:\n';
    errmsg = [errmsg,'   Batch input filenames must have the form "*_batch.m"'];
    errmsg = [errmsg,' where "*" begins with\n'];
    errmsg = [errmsg,'   "P"       for generation of Poisson problems\n'];
    errmsg = [errmsg,'   "CD"      for generation of convection-diffusion problems\n'];
    errmsg = [errmsg,'   "S"       for generation of Stokes problems\n'];
    errmsg = [errmsg,'   "NS"      for generation of Navier-Stokes problems\n'];
    errmsg = [errmsg,'   "itsolve" for iterative solution of linear systems.'];
    error('BATCH:err',errmsg);    
end  

% batch run
% switch to activate batch mode (off/on 0/1) (see "default.m")
BATCH=1;

% run appropriate driver
if strncmp(testproblem,'itsolve',7)
   load batchrun
   it_solve
   % save data
   gohome
   cd datafiles
   save batchrun_itsolve
else 
   if strncmp(testproblem,'P',1)  
      diff_testproblem
   elseif strncmp(testproblem,'CD',2)
      cd_testproblem
   elseif strncmp(testproblem,'S',1)
      stokes_testproblem   
   elseif strncmp(testproblem,'NS',2)
      navier_testproblem
   end
   % save data
   gohome
   cd datafiles
   save batchrun
end

% switch back to interactive mode
activemode
return
