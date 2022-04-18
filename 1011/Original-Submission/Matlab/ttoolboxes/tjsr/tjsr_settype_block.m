function [type]=tjsr_settype_block(type)
% This function belongs to tjsr
% Combines informations from type.block{:} and saves it unter type.
%
% Written by: tommsch, 2018

%compute JSR
val=cellfun(@(x) x.JSR, type.block, 'UniformOutput', 0); 
type.JSR=blockjsr(val{:});

%combine errortext
for i=1:type.counter.numblock
    %combine type.info.errortext
    if(isfield(type.block{i},'info')); 
        if(isfield(type.block{i}.info,'errortext') && ~isempty(type.block{i}.info.errortext));
            type.info.errortext=[type.info.errortext sprintf('\n') 'BLOCK:' num2str(i) sprintf('\n') type.block{i}.info.errortext]; %#ok<SPRINTFN>
        end
        if(isfield(type.block{i}.info,'infotext'));
            type.info.infotext=[type.info.infotext sprintf('\n') 'BLOCK:' num2str(i) sprintf('\n') type.block{i}.info.infotext]; %#ok<SPRINTFN>
        end
    end
end
    
%combine type.counter
type.counter.numstepsmall=0;
type.counter.numstepbig=0;
type.counter.iteration=0;
type.counter.totaltime=0;
type.counter.treetime=0;
type.counter.numberofvertex=0;
for i=1:type.counter.numblock
    if(isfield(type.block{i},'counter')); 
        if(isfield(type.block{i}.counter,'numstepsmall'));   type.counter.numstepsmall=   type.counter.numstepsmall+   type.block{i}.counter.numstepsmall;   end;
        if(isfield(type.block{i}.counter,'numstepbig'));     type.counter.numstepbig=     type.counter.numstepbig+     type.block{i}.counter.numstepbig;     end;
        if(isfield(type.block{i}.counter,'iteration'));      type.counter.iteration=      type.counter.iteration+      type.block{i}.counter.iteration;      end;
        if(isfield(type.block{i}.counter,'totaltime'));      type.counter.totaltime=      type.counter.totaltime+      type.block{i}.counter.totaltime;      end;
        if(isfield(type.block{i}.counter,'treetime'));       type.counter.treetime=       type.counter.treetime+       type.block{i}.counter.treetime;       end;
        if(isfield(type.block{i}.counter,'numberofvertex')); type.counter.numberofvertex= type.counter.numberofvertex+ type.block{i}.counter.numberofvertex; end;
    end
end
    
%combine type.cyclictree
type.cyclictree.ordering={};
type.cyclictree.smpflag=[];
type.cyclictree.v0={};
type.cyclictree.v0s={};
type.cyclictree.multiplicity=[];
type.cyclictree.orho=[];
type.cyclictree.maxlengthordering=[];
for i=1:type.counter.numblock
    if(isfield(type.block{i},'cyclictree')); 
        if(isfield(type.block{i}.cyclictree,'ordering'));          type.cyclictree.ordering=          [type.cyclictree.ordering          type.block{i}.cyclictree.ordering];          end;
        if(isfield(type.block{i}.cyclictree,'smpflag'));           type.cyclictree.smpflag=           [type.cyclictree.smpflag           type.block{i}.cyclictree.smpflag];           end;
        if(isfield(type.block{i}.cyclictree,'v0'));                type.cyclictree.v0=                [type.cyclictree.v0                type.block{i}.cyclictree.v0];                end;
        if(isfield(type.block{i}.cyclictree,'v0s'));               type.cyclictree.v0s=               [type.cyclictree.v0s               type.block{i}.cyclictree.v0s];               end;
        if(isfield(type.block{i}.cyclictree,'multiplicity'));      type.cyclictree.multiplicity=      [type.cyclictree.multiplicity      type.block{i}.cyclictree.multiplicity];      end;
        if(isfield(type.block{i}.cyclictree,'orho'));              type.cyclictree.orho=              [type.cyclictree.orho              type.block{i}.cyclictree.orho];              end;
        if(isfield(type.block{i}.cyclictree,'maxlengthordering')); type.cyclictree.maxlengthordering= [type.cyclictree.maxlengthordering type.block{i}.cyclictree.maxlengthordering]; end;
    end
end
    

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.