function type = tjsr_checkoptions(type);
% type = tjsr_checkoptions(type); 
% This function belongs to tjsr!
% make tests if options are sensible
% Input: 'type'-struct from tjsr
% Output: Messages on screen and in type
%
% Written by tommsch, 2018

% Start the computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Check if Gurobi is installed
    
    
    model.A = sparse([1 1 0; 0 1 1]);
    model.obj = [1 2 3];
    model.modelsense = 'Max';
    model.rhs = [1 1];
    model.sense = [ '<' '<'];
    params.outputflag=0;
    try 
        gurobiresult = gurobi(model,params);
        if(~isequal(gurobiresult.x,[ 1 0 1].')); 
            error('Gurobi seems to be installed, but does not work. Fix the Gurobi installation.');
        end        
    catch
        type.info.infotext=vprintf('Gurobi-solver not working. Fallback to matlab-linprog.\n','cpr','err','imp',[1,type.opt.verbose],'str',type.info.infotext);
    end
    
    if(type.opt.fastnorm>=2);
        type.info.infotext=vprintf('If fastnorm>=2, the intermediate bounds are wrong.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext);
        type.info.errortext=vprintf('If fastnorm>=2, the intermediate bounds are wrong.\n','str',type.info.errortext,'npr');        
    end
    if(type.opt.memory>=1);
        type.info.infotext=vprintf('Change parameters since ''memory'' is set.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext);
        type.info.errortext=vprintf('Change parameters since ''memory'' is set.\n','str',type.info.errortext,'npr');
        type.opt.plot='none';
    end
    if(type.opt.memory>=2)
        type.opt.testspectralradius=-inf;
        type.opt.testeigenplane=-inf;
    end
    if(type.opt.rholeqval && type.opt.fastnorm); 
        type.info.infotext=vprintf('You may like set <''fastnorm'',0> since ''rholeqval'' is set.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext);
        type.info.errortext=vprintf('You may like set <''fastnorm'',0> since ''rholeqval'' is set.\n','str',type.info.errortext,'npr');
        %vprintf('Press any key to continue.\n','cpr','err','imp',[0,type.opt.verbose]); pause;
    end;
    if(type.opt.alwaysout && type.opt.fastnorm); 
        vprintf('You should set <''fastnorm'',0> since ''alwaysout'' is set.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext); pause; 
    end;
    if(max(cellfun(@(x) tif(isempty(x),0,max(max(x))), type.cyclictree.ordering))>length(type.M_original)); 
        vprintf('Numbers given in ''ordering'' are larger than number of matrices.\n','cpr','err','imp',[0,type.opt.verbose]); 
        return; 
    end;
    if(type.opt.simplepolytope<0)
        type.info.infotext=vprintf('''simplepolytope'' should be greater-equal than zero, otherwise it will have no effect.\n','cpr','err','imp',[1,type.opt.verbose],'str',type.info.infotext);     
    end
    if(type.opt.testspectralradius>=0);        
        type.info.infotext=vprintf('''testspectralradius'' should be less than zero, to prevent false positives.\n','cpr','err','imp',[1,type.opt.verbose],'str',type.info.infotext); 
    end;
    if(type.opt.testspectralradius>=0);        
        type.info.infotext=vprintf('''testeigenplane'' should be less than zero, to prevent false positives.\n','cpr','err','imp',[1,type.opt.verbose],'str',type.info.infotext); 
    end;
%    if(type.opt.notestcycle && type.opt.delta>=1); 
%        type.info.infotext=vprintf('''notestcycle'' is set, but ''delta'' is geq 1. Set ''delta'' to 0.999.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext); type.opt.delta=0.999; 
%    end;
    if(type.opt.autoextravertex>=1); 
        type.info.infotext=vprintf('''autoextravertex'' should be smaller than 1.\n','cpr','err','imp',[1,type.opt.verbose],'str',type.info.infotext); 
    end;
    if(type.opt.alwaysout);      
        type.info.infotext=vprintf('!''alwaysout'' is set!\n','imp',[2,type.opt.verbose],'str',type.info.infotext); 
    end;
    if(type.opt.epsequal<0);     
        type.info.infotext=vprintf('epsequal<0. This will cause an error.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext); 
        type.info.errortext=vprintf('epsequal<0. This will cause an error.\n','str',type.info.errortext,'npr'); 
    end;
    %if(type.opt.epspolytope<type.opt.epslinprog);  
    %    type.info.infotext=vprintf('''epspolytope'' should be strictly greater than ''epslinprog''.\n','imp','cpr','err',[0,type.opt.verbose],'str',type.info.infotext); 
    %    type.info.errortext=vprintf('''epspolytope'' should be strictly greater than ''epslinprog''.\n','str',type.info.errortext,'npr'); 
    %end; 
    if(type.opt.epspolytope<0);  
        type.info.infotext=vprintf('epspolytope<0\n','imp',[1,type.opt.verbose],'str',type.info.infotext); 
    end; 
    if(type.opt.delta>1);        
        type.info.infotext=vprintf('delta>1\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext); 
        type.info.errortext=vprintf('delta>1.\n','str',type.info.errortext,'npr'); 
    end;
    if((~isempty(type.cyclictree.v0) || ~isempty(type.cyclictree.v0s) || ~isempty(type.cyclictree.extravertex) ||  ~isempty(type.cyclictree.smpflag)) && ~isequal(type.opt.invariantsubspace,'none'));
        type.info.infotext=vprintf('You should set <''invariantsubspace'',''none''> since you have given explicit starting vectors.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext); 
        type.info.errortext=vprintf('You should set <''invariantsubspace'',''none''> since you have given explicit starting vectors.\n','str',type.info.errortext,'npr'); 
    end; 
    if(type.opt.naturalselectiontype>=100 && type.opt.naturalselectiontype<1000);
        type.info.infotext=vprintf('''naturalselectiontype'' has a value where the algorithm will behave very badly.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext); 
        type.info.errortext=vprintf('''naturalselectiontype'' has a value where the algorithm will behave very badly.\n','str',type.info.errortext,'npr'); 
    end
    if(type.opt.naturalselection<0);
        type.info.infotext=vprintf('''naturalselection''<0, thus the algorithm cannot report intermediate bounds for the JSR.\n   If the algorithm still reports bounds, these are wrong!\n   If the algorithm reports an exact value, this value is correct.\n','cpr','err','imp',[0,type.opt.verbose],'str',type.info.infotext); 
        type.info.errortext=vprintf('''naturalselection''<0, thus the algorithm cannot report intermediate bounds for the JSR.\n   If the algorithm still reports bounds, these are wrong!\n   If the algorithm reports an exact value, this value is correct.\n','str',type.info.errortext,'npr'); 
    end
end