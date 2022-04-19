function varargout = makeorderinggraph(varargin)
% [ G ] = makeorderinggraph( oo, [options] )
% Constructs the graph given a partially ordered set.
%
% Input: 
%   oo      Matrix of orderings. Each column is one ordering OR
%           Cell array of orderings. Each cell consists of one column vector
%
% Options:
%   'value',val             Format of val: Nx2 cell array. First columns: name, second columns: column vector of values/cell array of strings
%   'fontsize',val          Fontsize of the vertex-labels
%   'labelonlyleaf'         Puts labels only on vertices without children
%   'labeldescription'      Puts the description (entries in 'value'{:,1}) 
%   'labelnumber'           Puts the number of the vertices
%   'labeledge',val         Labels the edges
%   'noplot'                No plot output. Just returns graph as object:
%   'verbose',val           Verbose level
%        
% Output:
%   G                       The graph as a matlab graph object
%   plot of the graph
%
% E.g.: oo=[1 1 1 2 1 2 2 2 1 1;0 2 1 0 2 1 2 1 1 2;0 0 0 0 1 0 0 1 2 1]; makeorderinggraph(oo,'labeledge');
%       findsmp(tgallery('rand_rho',2,2,10),'maxsmpdepth',10,'plot','tree','N',2)   
%       tjsr(tgallery('rand_rho',3,2,100),'plot','tree');
%       oo=randi(3,5,100)-1; makeorderinggraph(oo);
%
% Written by: tommsch, 2018 

% Changelog: 2019-03-30     added options 'labeldescription', 'labelnumber'
%                           cell arrays of text can be passed as vertex text

% XX Add cycles
% XX Add 'removedeadleafs'
% XX Add legend

    oo=varargin{1}; %size of oo changes during the computation
    varargin(1)=[];
    if(isnumeric(oo)); oo=num2cell(oo,1); end; %make oo to cell array
    
    [fontsize,varargin]=parsem('fontsize',varargin,5);
    [labelonlyleaf,varargin]=parsem('labelonlyleaf',varargin,0);
    [labeledge,varargin]=parsem('labeledge',varargin,0);
    [labeldescription,varargin]=parsem('labeldescription',varargin);
    [labelnumber,varargin]=parsem('labelnumber',varargin);
    [verbose,varargin]=parsem({'verbose','v'},varargin,1);    
    [value,varargin]=parsem('value',varargin,{}); %format: 'name',val
    [noplot,varargin]=parsem('noplot',varargin);
    
    parsem(varargin,'test');
    
    
%Make input checks    
    nvalue=size(value,1);
    if(nvalue>0)
        if(~iscell(value)); error('Argument of ''value'' must be a Nx2 cell'); end;
        if(~isequal(size(value,2),2)); error('Argument of ''value'' must be a Nx2 cell'); end;

        for i=1:nvalue
            if(~ischar(value{i,1})); error('First entries of value must be a string.'); end;
            if(isrow(value{i,2})); value{i,2}=value{i,2}.'; end;
            if(~isequal(size(value{i,2}),[size(oo,2),1])); error('Second entries must be a column vector with as many values as elements in oo.'); end;
        end
    end

%Preprocess
    vprintf('Construct Graph with %i vertices. ',numel(oo),'imp',[1,verbose]);
    
    G = zeros(3,numel(oo)); %the nodes (vertices)
    
    %if there is no empty ordering, add it
    if(~isempty(oo) && any(oo{1})); 
        oo=horzcat(1,oo); oo{1}=[]; %add empty ordering (gets the number zero in the next steps)
        for i=1:nvalue; value{i,2}=[0; value{i,2}]; end; %add values for empty ordering
    end
    noo=numel(oo);
    
    for i=1:noo;
        if(isrow(oo{i})); oo{i}=oo{i}.'; end; %make to column vector
        oo{i}=[0; removezero(oo{i},1)]; %add the root vertex temporarily
    end;
    
    for i=1:noo;  %construct graph
        [found,idx]=searchincellarray(oo{i}(1:end-1),oo,1);
        if(found);
            G(:,i)=[idx(1); i; oo{i}(end)]; %search for vertices: Format: node, node, oo-idx
        end
    end
    G=removezero(G,2); %remove zeros in between
        
%Pass object to matlab and let matlab do some work
    vprintf('Pass object to matlab. ','imp',[1,verbose]);
    if(isempty(G));
        if(nargout==1); varargout{1}=G; end;
        return;
    else
        G=digraph(G(1,:),G(2,:),G(3,:)); 
        for i=1:nvalue; eval(['G.Nodes.' value{i,1} '=value{i,2};']); end; %Add values to the table
    end

%Plot Graph    
    vprintf('Plot Graph. ','imp',[1,verbose]);
    if(~noplot)
        h=plot(G,'Layout','layered','NodeLabel',[],'AssignLayers','asap','Direction','r');
    end
    
%Add text    
    vprintf('Add Text. ','imp',[1,verbose]);
    vertextext=cell(size(G.Nodes,1),1); 
    for i=1:length(h.XData); vertextext{i}=''; end; %make empty text
    for i=1:length(h.XData); 
        for j=1:nvalue
            val=value{j,2}(i);
            if(iscell(val) && ischar(val{1})); val=val{1}; %do nothing
            elseif(isnumeric(val)); val=num2str(val.');
            end;
            if(labeldescription); val=[value{j,1} ': ' val]; %#ok<AGROW>
            end
            vertextext{i}=[vertextext{i} val newline]; 
        end; 
        if(labelnumber); vertextext{i}=[num2str(i) ': ' newline vertextext{i}];
        end
    end
    for i=1:length(h.XData); 
        if(labelonlyleaf && any(G.Edges.EndNodes(:,1)==i)); 
            continue; end;
        text(h.XData(i)+0.1,h.YData(i),vertextext{i},'FONTSIZE',fontsize); 
    end; 

%Process Edges
    vprintf('Process Edges. ','imp',[1,verbose]);
    J=max(cellfun(@max,oo)); %number of different matrices in tree
    LINESTYLE={'-',':','--','-.'};
    for i=1:J;
        idx=G.Edges.Weight==i;
        if(labeledge); 
            labeledge(h,G.Edges.EndNodes(idx,1),G.Edges.EndNodes(idx,2),num2str(i)); 
        end;
        highlight(h,G.Edges.EndNodes(idx,1),G.Edges.EndNodes(idx,2),'EdgeColor',num2color(i),'LineWidth',1.5,'LineStyle',LINESTYLE{mod(i,4)+1});
    end;

%Postprocessing
    if(~isempty(G.Nodes)); 
        highlight(h,1,'NodeColor','k','MarkerSize',12); %Highlight root
    end;
    axis off
    if(nargout==1); varargout{1}=G; end;
    vprintf('\n','imp',[1,verbose]);

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 