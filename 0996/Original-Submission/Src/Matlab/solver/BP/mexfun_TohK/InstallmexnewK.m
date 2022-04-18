%%
%% Run this script in Matlab command window 
%%

   function InstallmexnewK(recompile)

   if (nargin==0); recompile = 0; end

   curdir = pwd;  
   fprintf(' current directory is:  %s\n',curdir);    
%%
%% generate mex files in Mexfun
%% 
   computer_model = computer;
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
   tmp = version('-release'); 
   matlabrelease = str2num(tmp(1:4));
%%
   if strcmp(computer_model,'PCWIN')
      str0 = ['''',matlabroot,'\extern\lib\win32\lcc\'' '];
      if (exist(eval(str0),'dir')==7)
         str1 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwblas.lib''  '];
      else
         str1 = ['''',matlabroot,'\extern\lib\win32\microsoft\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win32\microsoft\libmwblas.lib''  '];
      end
      libstr = [str1,str2];
   elseif strcmp(computer_model,'PCWIN64')
      str0 = ['''',matlabroot,'\extern\lib\win64\lcc\'' '];
      if (exist(eval(str0),'dir')==7)
         str1 = ['''',matlabroot,'\extern\lib\win64\lcc\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win64\lcc\libmwblas.lib''  '];
      else
         str1 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwblas.lib''  '];
      end
      libstr = [str1,str2];
   else
      libstr = '  -lmwlapack -lmwblas  '; 
   end
   mexcmd = 'mex -O  -largeArrayDims  -output ';    
%%
   if (matlabversion < 7.3)  && (matlabrelease <= 2008)
      error(' needs MATLAB version 7.4 and above'); 
   end
   fsp = filesep;   
%%
%%
%%
   src = curdir;
   fprintf ('\n Now compiling the mexFunctions in:\n'); 
   fprintf (' %s\n',src);    

%    fname{1} = 'mexbwsolveK'; 
%    fname{2} = 'mexfwsolveK';
%    fname{3} = 'mexFnormK'; 
%    fname{4} = 'mexmatK';
%    fname{5} = 'mexnnzK';
%    fname{6} = 'mexprojsocK';
%    fname{7} = 'mexsmatK'; 
%    fname{8} = 'mexsvecK';
%    fname{9} = 'mextriangK';
%    fname{10}= 'mexeigPartialK';
%    fname{11}= 'mexeigK'; 
   
   fname{1} = 'mexFnormK'; 
   fname{2}= 'mexeigPartialK';
   fname{3}= 'mexeigK'; 

%%
   ext = mexext; 
   for k = 1:length(fname)
      existmex = exist([fname{k},'.',ext]); 
      if (existmex ~= 3) || (recompile)
         cmd([mexcmd,fname{k},'  ',fname{k},'.c  ',libstr]);  
      end
   end 
   
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************
