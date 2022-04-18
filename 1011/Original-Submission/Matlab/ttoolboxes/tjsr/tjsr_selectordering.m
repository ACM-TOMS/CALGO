function [v0, v0s, oo, smpflag, mult, type] = tjsr_selectordering(v0, v0s, oo, smpflag, mult, type);
    %checks if there are multiple eigenvectors, and put same eigenvectors into the same class
    %repeats orderings and cycles orderings to find more matches
    %returns everything as row vector   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%These are values for a non-trivial example to test the algorithm
    % Ex1:
    % % tjsr({[2 1; 0 -2],[2 0; -1 -2]},'invariantsubspace','none','maxsmpdepth',1,'verbose',2)
    % % tjsr(tgallery('ex_rand_rho',2,2,7),'invariantsubspace','none','verbose',3,'nearlycandidate',0.99,'maxsmpdepth',10,'nobalancing')
    %
    % Ex2:
%     v0_mult={}; nv0_mult={}; o_mult={}; no_mult={};  %DEBUG
%     v0_mult{1}=[1;1];v0_mult{10}=[1;1];
%     v0_mult{2}=[2;2];v0_mult{9}=[2;2];v0_mult{6}=[2;2];
%     v0_mult{13}=[3;3];v0_mult{4}=[3;3];v0_mult{7}=[3;3];
%     v0_mult{8}=[4;4];v0_mult{11}=[4;4];
%     v0_mult{12}=[5;5];
%     v0_mult{5}=[6;6];
%     v0_mult{3}=[7;7];
%     nv0_mult{1}=[-1;-1];
%      
%     o_mult{1}=[2 1]; o_mult{2}=[1 2]; o_mult{13}=[1 2 1 2];
%     o_mult{4}=[1 3 3]; o_mult{5}=[3 3 1]; o_mult{6}=[3 1 3];
%     o_mult{7}=[2]; o_mult{8}=[2 2];
%     o_mult{9}=[1 3]; o_mult{10}=[1 3 1 3]; o_mult{11}=[3 1 3 1]; o_mult{12}=[3 1];
%     o_mult{3}=[3];
%     no_mult{1}=[1 2 3];
    %
    % Ex3:
    
     %v0={};  oo={};   %DEBUG
     %oo{1}=[1 1 1 1 2]'; oo{2}=[1 1 1 2 1]'; oo{3}=[1 1 2 1 1]'; oo{4}=[1 2 1 1 1]'; oo{5}=[2 1 1 1 1]'; 
     %oo{6}=[1 2]'; oo{7}=[1 2]';
     %oo{8}=[1]'; oo{9}=[2]'; oo{10}=[2]';
     %v0{1}=[1;1]; v0{2}=[2;2];v0{3}=[3;3]; v0{4}=[4;4];v0{5}=[5;5];
     %v0{6}=[1;1]; v0{7}=[2;2];
     %v0{8}=[3;3]; v0{9}=[6;6]; v0{10}=[10;10];
     %smpflag=[0 0 0 0 0 1 1 1 1 2];
     %mult=[1 1 1 1 1 2 2 3 3 3];
     %type.opt.verbose=3; type.info.infotext=''; type.opt.epsequal=eps; type.opt.noclassify=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(any(cellfun(@(x) size(x,2),oo)>1)); 
        return; end;
    if(size(oo,2)==1); 
        return; end;
    o_short = reducelength(oo);
    
    type.info.infotext = vprintf('   Compute classes of eigenvectors and orderings (clusters) - ','imp',[3,type.opt.verbose],'str',type.info.infotext);
    val = [v0{:}]; 
    val = [real(val); imag(val)]; 
    val = val+2*abs(min(min(val)))+1; %workaround for strange behaviour of uniquetol: if one entry is zero, uniquetol behaves like unique
    [~,idx_vs_v0class,idx] = uniquetol( val', type.opt.epsequal, 'ByRows',true );
    [~,idx_vs_v0class] = sort( idx_vs_v0class ); 
    [~,idx_vs_v0class] = sort( idx_vs_v0class ); %remove "holes" in the numbering off class
    idx_vs_v0class = idx_vs_v0class(idx).';      %if v0class{i}==v0class{j} then v0{i}==v0{j}

    [~,idx_vs_o,idx] = uniquecell(oo);
    [~,idx_vs_o] = sort(idx_vs_o); 
    [~,idx_vs_o] = sort(idx_vs_o); %remove "holes" in the numbering off class
    idx_vs_o = idx_vs_o(idx).'; %if oclass{i}==oclass{j} then allo(i)==allo(j)
    
    [~,idx_vs_oclass,idx] = uniquecell( o_short );
    [~,idx_vs_oclass] = sort( idx_vs_oclass ); 
    [~,idx_vs_oclass] = sort( idx_vs_oclass ); %remove "holes" in the numbering off class
    idx_vs_oclass = idx_vs_oclass( idx ).'; %if oclass{i}==oclass{j} then reducelength(i)==reducelength(j)
        
    type.info.infotext=vprintf('   %i eigenvector(s), %i ordering(s) found.\n',max(idx_vs_v0class), max(idx_vs_oclass),'imp',[3,type.opt.verbose],'str',type.info.infotext);
    
    %partitionate problem in subproblems
    type.info.infotext = vprintf('   Compute connectivity between clusters, partitionate problem in subproblems - ','imp',[3,type.opt.verbose],'str',type.info.infotext);
    v0class_vs_oclass = cell(1,max(idx_vs_v0class)); 
    for i = 1:length( v0class_vs_oclass ); 
        v0class_vs_oclass{i} = sort( idx_vs_oclass(idx_vs_v0class==i) ); 
        v0class_vs_oclass{i} = unique( v0class_vs_oclass{i} ); end; 
    v0class_vs_o = cell( 1, max(idx_vs_v0class) ); 
    for i = 1:length( v0class_vs_o ); 
        v0class_vs_o{i} = sort( idx_vs_o(idx_vs_v0class==i) ); 
        v0class_vs_o{i} = unique( v0class_vs_o{i} ); end; 
    
    
    %look which v0-class hits which ordering-class    
    v0class_vs_oclass_2=v0class_vs_oclass; %oclass_2 is actually the same as "problem"
    for i=1:length(v0class_vs_oclass_2); 
        for j=i+1:length(v0class_vs_oclass_2) %#ok<ALIGN>
            if(intersect(v0class_vs_oclass_2{i},v0class_vs_oclass_2{j}));
                v0class_vs_oclass_2{i}=union(v0class_vs_oclass_2{i},v0class_vs_oclass_2{j}); 
                v0class_vs_oclass_2{j}=v0class_vs_oclass_2{i}; end; %union all v0_classes
    end; end;

    %each different v0_class_2 is a problem
    [~,~,v0class_vs_problem]=uniquecell(v0class_vs_oclass_2); 
    v0class_vs_problem=v0class_vs_problem';
    type.info.infotext=vprintf('   %i subproblem(s) found.\n',max(v0class_vs_problem),'imp',[3,type.opt.verbose],'str',type.info.infotext);
    
    problem_vs_oclass=cell(1,max(v0class_vs_problem)); 
    problem_vs_v0class=cell(1,max(v0class_vs_problem));
    for i=1:max(v0class_vs_problem); 
        problem_vs_oclass{i}=unique(cell2mat(v0class_vs_oclass(v0class_vs_problem==i)));    %compute which ordering-indices belong to the subproblems
        problem_vs_v0class{i}=find(v0class_vs_problem==i); end;                         %compute which v0-indices belong to the subproblems
    
    if(type.opt.noclassify)
        idx_all=num2cell(1:length(v0));
    else
        %solve problem for each sub-problem
        %XX if there are too many candidates, we cant compute all combinations and thus need to make some trivial choice
        type.info.infotext=vprintf('   Choose candidates for each subproblem.\n','imp',[3,type.opt.verbose],'str',type.info.infotext);
        idx_v0class=cell(size(problem_vs_oclass));
        for i=1:max(v0class_vs_problem);
            %make all possible choices of eigenvectors
            v0choice={}; 
            for j=1:numel(problem_vs_oclass{i});
                v0choice=[v0choice num2cell(combnk(problem_vs_v0class{i},j),2)']; end; %#ok<AGROW>
            type.info.infotext=vprintf('   Size of subproblem no. %i: %i.',i,length(v0choice),'imp',[3,type.opt.verbose],'str',type.info.infotext);

            %test if choice covers whole problem
            ocovered=cell(1,length(v0choice));
            for j=1:length(v0choice); 
                for k=1:length(v0choice{j}); 
                    ocovered{j}=[ocovered{j}, idx_vs_oclass(idx_vs_v0class(:)==v0choice{j}(k))]; end; %compute which clusters are covered
                ocovered{j}=unique(ocovered{j}); end; 
            idx=find(cellfun(@(x) isequal(x,problem_vs_oclass{i}), ocovered)); %save which choices cover everything

            %choose the choice with minimal length  of orderings involved
            type.info.infotext=vprintf('   Compute weights.','imp',[3,type.opt.verbose],'str',type.info.infotext);        
            CUMWEIGHT=inf*ones(size(v0choice)); %initialize variable, set to infinity so that I can use min(WEIGHT) afterwards
            for j=idx %go through every possible choice                
                CUMWEIGHT(j)=0; %set to zero, since I want to add up numbers
                for k=1:length(v0choice{j});
                    CAND=idx_vs_v0class==v0choice{j}(k); %choose candidates
                    CUMWEIGHT(j)=CUMWEIGHT(j)+min(cellfun(@length, oo(CAND))); end; %add length of candidates to the WEIGHT
            end;
            idx=intersect(idx,find(CUMWEIGHT==min(CUMWEIGHT)));

            %choose the idx with minimal ordering-length of v0choice(idx)
            LENGTH=(cellfun(@length, v0choice(idx)));
            idx=idx(LENGTH==min(LENGTH));

            %choose the idx with least ordering-weight of v0choice(idx) %XX Test what happens with multiple leading ev
            WEIGHT=cellfun(@numel ,v0choice(idx));
            idx=idx(WEIGHT==min(WEIGHT));
            allv0_idx=[v0choice{idx(1)}]; %make idx point to (original) allv0

            idx_v0class{i}=allv0_idx(:).'; %choose all minimal lengths
            type.info.infotext=vprintf('   Subproblem solved.\n','imp',[3,type.opt.verbose],'str',type.info.infotext);
        end;
        idx_v0class=[idx_v0class{:}];
        idx_all=arrayfun(@(x) find(idx_vs_v0class==x),idx_v0class,'UniformOutput',0);
    end
    
    
    v0=cellfun(@(x) v0(x(1)),idx_all);
    v0s=cellfun(@(x) v0s(x(1)),idx_all);
    
    val=oo;
    %put all orderings together
    for i=1:length(idx_all)
        val2=uniquecell(reducelength(val(idx_all{i}),1));
        val2=removecombination(val2);
        maxSize = max(cellfun(@numel, val2));                   % Get the maximum vector size
        padd = @(x) [x; zeros(maxSize-numel(x),1)];             % Create an anonymous function
        oo{i} = cellfun(padd, val2, 'UniformOutput', false);    % Pad each cell with 0
        oo{i} = cell2mat(oo{i}); end;                           % Vertically concatenate cells
    
    oo(length(idx_all)+1:end)=[];
    
    mult=cellfun(@(x) max(mult(x)),idx_all);
    smpflag=cellfun(@(x) min(smpflag(x)),idx_all);

    type.info.infotext=vprintf('   Everything done.\n','imp',[3,type.opt.verbose],'str',type.info.infotext);
end