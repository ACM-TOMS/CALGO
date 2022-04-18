function tjsr_printerrortext(type,alltype)
% This function belongs to tjsr
% Prints out information in the fields .errortext
%
% Written by: tommsch, 2019

% XX Function is quite untested and should be tested

try
    if( type.opt.verbose>=0 && ~isempty(type.info.errortext));
        if(~isempty(alltype)); 
            for i=1:length(alltype); 
                if(~isempty(alltype{i}.info.errortext));
                    vprintf('\nImportant Messages (Run: %i) \n================\n%v\n',i,alltype{i}.info.errortext,'cpr','err');
                end
            end
        else
            if(type.counter.numblock>1)
                for i=1:type.counter.numblock
                    if(isfield(type.block{i},'info') && isfield(type.block{i}.info,'errortext') && ~isempty(type.block{i}.info.errortext))
                        vprintf('\nImportant Messages (Block: %i) \n================\n%v\n',i,type.info.errortext,'cpr','err');
                    end
                end
            else
                if(~isempty(type.info.errortext));
                    vprintf('\nImportant Messages: \n================\n%v\n',type.info.errortext,'cpr','err');
                end
            end
        end
    end 
catch
    vprintf('\nError while printing error-messages from tjsr.\n','cpr','err');            
end

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.