function type = tjsr_saveoutput(type,savelevel)
% type = tjsr_saveoutput(type,savelevel)
% Saves most of the output of the running algorithm in files
% This function belongs to tjsr!
%
% Input:
%       type
%       savelevel       how much shall be saved
%                           >=3 a lot
%                           <3  not much
%
% Written by tommsch, 2018

%#ok<*TRYNC>
    type.info.infotext=vprintf('Write to diary.\n','imp',[3,type.opt.verbose],'str',type.info.infotext); 
    try
        if(~isequal(0,type.opt.diary)); 
           diary; diary;
        end;
    end

    try
        if(type.opt.plot); 
            print(['plot tjsr ' datestr(datetime)],'-dpng')
        end
    end
    try 
        if(savelevel >= 3 && type.opt.profile); 
            type.info.infotext=vprintf('Write profiling-info to diary. This may take a while.\n','imp',[2,type.opt.verbose],'str',type.info.infotext); 
            profsave(profile('info'),['profile tjsr ' datestr(datetime)]); 
            profile on;
        end;    
    end
    try
        if(savelevel>=3)
            type.info.infotext=vprintf('Write ''type'' to diary.\n','imp',[3,type.opt.verbose],'str',type.info.infotext); 
            save(['type tjsr ' datestr(datetime)],'type');
        end
    end

    
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 