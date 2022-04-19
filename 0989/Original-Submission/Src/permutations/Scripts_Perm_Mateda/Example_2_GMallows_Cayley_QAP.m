% This script contains an example of execution of EDAs for solving permutation problems using the toolbox presented in the paper entitled:
% perm_mateda: A matlab toolbox of estimation of distribution algorithms for permutation-based combinatorial optimization problems
% submitted to the journal ACM Transactions on Mathematical Software.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- Using Generalized Mallows-Kendall EDA for solving the Permutation Flowshop Scheduling Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 global QAPInstance                   % The selected QAP instance is saved in the global variable QAPInstance
 
 global FitnessImprovement
 
  ReadQAPInstance('../Problems/Instances/QAP/tai40b');       % QAP (Taillard Benchmark -> tai40b.dat, n=40)a

 [distance, flow, n] = QAPInstance{:};                           % QAP     Problem
 
 nrepetitions=30; 	% Number of repetitions of the algorithm
 ngen = 500;                   % Number of generations 
 NumbVar = n;                  % Number of variables
 PopSize = 10*NumbVar;         % Population size
 F = 'EvalQAP';               % Optimization problem     [It should be modified according to the problem]
 cache  = [0,0,1,0,0];         % Array to indicate which information
 
 Card = [ones(1,NumbVar);NumbVar*ones(1,NumbVar)];            % All variables have values between 1 and NumbVar
 edaparams{1} = {'seeding_pop_method','InitPermutations',{}}; % The initial population is uniquely composed of permutations 
 edaparams{2} = {'learning_method','GMallows_cayley_learning',{0.001,10,100,'setMedianPermutation',0.1}};          % Parameters of the learning algorithm
 edaparams{3} = {'sampling_method','GMallows_cayley_sampling',{PopSize-1,1}};                                    % Parameters of the sampling algorithm
 edaparams{4} = {'replacement_method','pop_aggregation_theta',{'fitness_ordering'}};                                    % Pop. aggregation is used as replacement method, Pop + sampled pop 
 selparams(1:2) = {NumbVar/PopSize,'fitness_ordering'};                                                          % Parameters used for selection (Truncation parameter and max_gen)
 edaparams{5} = {'selection_method','truncation_selection',selparams};                                           % Selection method used
 edaparams{6} = {'stop_cond_method','max_gen',{ngen}};                                                           % The algorithm stops when the max number of generations is reached
 
%%%% Run EDA executions and collect optimization statistics, execution times, and the models learn at each generation.
 
 Expe_AllStat = [];
 Expe_AllCache = [];
 Expe_AllTimes = {}
 
 for j=1:nrepetitions,
   FitnessImprovement = 0;
   [AllStat,Cache]=RunEDA(PopSize,NumbVar,F,Card,cache,edaparams);   
   for i=1:ngen, 
     Expe_AllStat(j,i) = AllStat{i,1}(2);
     Expe_AllTimes{j}(i,:) = AllStat{i,6};
     Expe_AllCache(j,i) = mean(Cache{3,i}{3});
   end
 end
 
 %%% Save results in file.
 fname = ['Results_GMallows_Cayley_QAP.mat']
 eval(['save ',fname, ' Expe_AllCache Expe_AllStat Expe_AllTimes'])
 
