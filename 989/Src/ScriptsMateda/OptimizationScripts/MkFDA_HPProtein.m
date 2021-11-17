% EXAMPLE 6:  Markov Chain FDA for the HP protein model.

 global HPInitConf;   % This is the HP protein instance, defined as a sequence of zeros and ones
 HPInitConf =  [zeros(1,12),1,0,1,0,1,1,0,0,1,1,0,0,1,1,0,1,1,0,0,1,1,0,0,1,1,0,1,1,0,0,1,1,0,0,1,1,0,1,0,1,zeros(1,12)]; 
 % The number of variables is equal to the sequence length and each
 % variables takes values in {0,1,2}
 PopSize = 1000; NumbVar = size(HPInitConf,2); cache  = [1,1,1,1,1]; Card = 3*ones(1,NumbVar);   maxgen = 10;
 % The Markov chain model(Cliques) is constructed specifying the number of
 % conditioned (previous) variables. In the example below this number is
 % 1., i.e. p(x) = p(x0)p(x1|x0) ... p(xn|xn-1) 
 Cliques = CreateMarkovModel(NumbVar,2);
 
 F = 'EvaluateEnergy'; % HP protein evaluation function
 edaparams{1} = {'learning_method','LearnFDA',{Cliques}};
 edaparams{2} = {'sampling_method','SampleFDA',{PopSize}};
 edaparams{3} = {'repairing_method','HP_repairing',{}}; % Repairing method used to guarantee that
                                                        % solutions do not self-intersect
 edaparams{4} = {'stop_cond_method','max_gen',{maxgen}};
 [AllStat,Cache]=RunEDA(PopSize,NumbVar,F,Card,cache,edaparams) 
 
 % To draw the resulting solution use function PrintProtein(vector),
 % where vector is the best solution found.
 vector = AllStat{maxgen,2}
 PrintProtein(vector)