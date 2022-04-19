% EXAMPLE 15:  Continuous EDAs that learn mixtures of distributions
 %              for  the trajectory problem (see previous examples for
 %              details on this problem)

 global MGADSMproblem
 cd ../trajectory % The spacecraft-trajectory problem instance is in this directory
 load EdEdJ
 cd ../Mateda2.0
 
 NumbVar = 12;
 PopSize = 5000; 
 F = 'EvalSaga';
 Card(1,:) = [7000,0,0,0,50,300,0.01,0.01,1.05,8,-1*pi,-1*pi];
 Card(2,:) = [9100,7,1,1,2000,2000,0.90,0.90,7.00,500,pi,pi]; 
 cache  = [0,0,1,0,1]; 
  
 learning_params(1:5) = {'vars','ClusterPointsKmeans',10,'sqEuclidean',1};
 edaparams{1} = {'learning_method','LearnMixtureofFullGaussianModels',learning_params};
 edaparams{2} = {'sampling_method','SampleMixtureofFullGaussianModels',{PopSize,3}};
 edaparams{3} = {'replacement_method','best_elitism',{'fitness_ordering'}};
 selparams(1:2) = {0.1,'fitness_ordering'};
 edaparams{4} = {'selection_method','truncation_selection',selparams};
 edaparams{5} = {'repairing_method','SetWithinBounds_repairing',{}};
 edaparams{6} = {'stop_cond_method','max_gen',{5000}};
 [AllStat,Cache]=RunEDA(PopSize,NumbVar,F,Card,cache,edaparams) 
