
% This scripts compares crossover and mutation algorithms in terms
% of the average fitness reached in each generation
% First, nruns runs of each EDAs are executed and the statistics of the
% algorithm are saved in a file ResultsCXMethod_i_Exp_j.mat
% where i is one of the crossover methods and j goes from 1 to nruns.
% In a second step, the files are read and analyzed to create the curves
% that describe the behavior of the algorithms.


PopSize = 300; 
n = 40; % Number of variables
cache  = [1,1,1,1,1]; % Save all the information about the models
Card = 2*ones(1,n);   % All variables are binary
maxgen = 50;          % Maximum number of generations
F = 'sum';% Onemax function;

% In order to achieve the desired factorization we use a junction tree. With 0 as parameter we have not overlappings. 
% Therefore, the cliques have a unique variable and they are independent.
Cliques = CreateMarkovModel(n, 0);

edaparams{1} = {'selection_method','truncation_selection',{0.5,'fitness_ordering'}};
edaparams{2} = {'stop_cond_method','max_gen',{maxgen}};

for i=2:2, % Modify to include all (the 5) selection operators 
 if i==1  
  edaparams{3} = {'learning_method','LearnFDA',{Cliques}};      % UMDA
  edaparams{4} = {'sampling_method','SampleFDA',{PopSize}};
 elseif i==2
  edaparams{3} = {'learning_method','LearnTwoPointCrossover',{PopSize}};      % Two point crossover
  edaparams{4} = {'sampling_method','SampleTwoPointCrossover',{PopSize}}; 
 end

 for Exp=[1:5],     
    filename = ['ResultsCXMethod_',num2str(i),'_Exp_',num2str(Exp),'.mat'];     
    [AllStat,Cache]=RunEDA(PopSize,n,F,Card,cache,edaparams); 
    eval(['save ', filename,  ' AllStat Cache']);
 end
end
  

% Extraction of the information from the files

for i=1:5, % Modify to include all selection operators   
  for Exp=[1:5],     
    filename = ['ResultsCXMethod_',num2str(i),'_Exp_',num2str(Exp),'.mat'];     
    eval(['load ', filename]);  % The file for selection operator i and experiment Exp is loaded
    for j=1:maxgen,
     MaxFitness{i}(Exp,j) = AllStat{j,1}(1);  % Maximum fitness of the population at each generation for each experiment
     MeanFitness{i}(Exp,j) = AllStat{j,1}(2); % Mean fitness of the population at each generation for each experiment    
    end,
  end
end

% Creation of two figures showing the Maximum and Mean fitness of the
% population in each generation

thecolors='rcbgk';
theforms ='+*ov<';
close all
fs = 14;

h=figure
for i=[1:5] % Modify to include all crossover and/or mutation operators  
  plot([1:2:maxgen],mean(MaxFitness{i}(:,1:2:maxgen)),[thecolors(i),theforms(i)],'LineWidth',2);
  hold on
end
G = ylabel('Maximum fitness'); 
H = xlabel('Generations');
legend(' truncation  ','proportional',' Tournament ',' LinearRank ','NonLinearRank','Location','NorthEast');
set(G,'Fontsize',fs);
set(H,'Fontsize',fs);
filename = ['MaxfitnessCX.eps'];
saveas(h,filename,'psc2');

h=figure
for i=[1:5] % Modify to include all crossover and/or mutation operators  
  plot([1:2:maxgen],mean(MeanFitness{i}(:,1:2:maxgen)),[thecolors(i),theforms(i)],'LineWidth',2);
  hold on
end
G = ylabel('Mean fitness'); 
H = xlabel('Generations');
legend(' truncation  ','proportional',' Tournament ',' LinearRank ','NonLinearRank','Location','NorthEast');
set(G,'Fontsize',fs);
set(H,'Fontsize',fs);
filename = ['MeanfitnessCX.eps'];
saveas(h,filename,'psc2');
 