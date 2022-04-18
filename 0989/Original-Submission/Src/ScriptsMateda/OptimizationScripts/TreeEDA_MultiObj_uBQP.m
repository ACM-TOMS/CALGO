%% EXAMPLE : MULTIO-OBJECTIVE binary uBQP
%Tree-EDA  for multi-objective uBQP problem

PopSize = 300; n = 100;
cache  = [0,0,1,1,1];
Card = 2*ones(1, n);
maxgen = 15;

global all_obj
[all_obj] = LoadUBQInstance('HardInst1__n_100.dat',100,2)
%[all_obj] = LoadUBQInstance('Art_type_3_n_100.dat',100,2)

F = 'Eval_uBQP'; % uBQP function
%F = 'Eval_MultTrap'; % Multiple deceptive function
%F = 'Eval_MultTrap_2'; % Multiple deceptive function
selparams(1:2) = {0.5,'ParetoRank_ordering'};

edaparams{1} = {'stop_cond_method','max_gen',{maxgen}};
edaparams{2} = {'replacement_method','best_elitism',{'ParetoRank_ordering'}};
edaparams{3} = {'learning_method','LearnTreeModel',{}};
edaparams{4} = {'sampling_method','SampleFDA',{PopSize}}; 

% Launch EDA

nexp = 30;
 AllMat = zeros(n);
 for i=2:nexp,
   [AllStat,Cache]=RunEDA(PopSize,n,F,Card,cache,edaparams)
   mat = zeros(n);  
   for j=1:size(Cache,2)     
     for k=1:n,
        if Cache{3,j}{1}(k,4)>0 
          mat(Cache{3,j}{1}(k,3),Cache{3,j}{1}(k,4)) = 1;  
          mat(Cache{3,j}{1}(k,4),Cache{3,j}{1}(k,3)) = 1;
        end  
     end
     AllMat = AllMat + mat;
   end   
 end  





