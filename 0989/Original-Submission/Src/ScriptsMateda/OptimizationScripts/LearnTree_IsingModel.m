 % EXAMPLE :  
 % Tree EDA for the Ising Model using MPE-sampling (i.e.
 % the most probable individual from the net is added to the population
 % during sampling). The algorithm stops when the optimum is found or the
 % max number of generations is reached
 
 
 global Isinglattice;
 global Isinginter;
 PopSize = 500; n = 64; cache  = [0,0,0,0,0]; Card = 2*ones(1,n); 
 F = 'EvalIsing'; MaxGen = 150;  MaxVal = 128;
 [Isinglattice, Isinginter] = LoadIsing(n, 6);
 stop_cond_params = {MaxGen,MaxVal};
 edaparams{1} = {'stop_cond_method','maxgen_maxval',stop_cond_params};
 edaparams{2} = {'learning_method','LearnTreeModel',{}};
 edaparams{3} = {'sampling_method','SampleFDA',{PopSize}}; 
 %edaparams{2} = {'sampling_method',' SampleMPE_BN',{PopSize}};
 [AllStat,Cache]=RunEDA(PopSize,n,F,Card,cache,edaparams);
   
 
 
 
 PopSize = 200; n = 20; cache  = [1,1,1,1,1]; Card = 2*ones(1,n); 
 F = 'ZeroOrOne'; MaxGen = 150; MaxVal = n;

 stop_cond_params = {MaxGen,MaxVal};
 edaparams{1} = {'stop_cond_method','maxgen_maxval',stop_cond_params};
 edaparams{2} = {'learning_method','LearnTreeModel',{}};
 edaparams{3} = {'sampling_method','SampleFDA',{PopSize}}; 
 %edaparams{2} = {'sampling_method',' SampleMPE_BN',{PopSize}};
 tt = zeros(1,20);
 aa = zeros(1,20);
 dd = zeros(1,20);
 for j=1:100,
  [AllStat,Cache]=RunEDA(PopSize,n,F,Card,cache,edaparams);
  dd(j) = sum(AllStat{end,2}); %mean(mean(Cache{1,1}) - mean(Cache{2,1}));
  tt = tt + mean(Cache{2,1});
  aa = aa + mean(Cache{1,1});
 end
 plot(mean(Cache{2,1}),'.')
 
 
 
  PopSize = 200; n = 20; cache  = [1,1,1,1,1]; Card = 2*ones(1,n); 
 F = 'ZeroOrOne'; MaxGen = 150; MaxVal = n;
 Cliques = CreateMarkovModel(n,2);
 stop_cond_params = {MaxGen,MaxVal};
 edaparams{1} = {'stop_cond_method','maxgen_maxval',stop_cond_params};
 edaparams{2} = {'learning_method','LearnFDA',{Cliques}};
 edaparams{3} = {'sampling_method','SampleFDA',{PopSize}};
 [AllStat,Cache]=RunEDA(PopSize,n,F,Card,cache,edaparams);
 