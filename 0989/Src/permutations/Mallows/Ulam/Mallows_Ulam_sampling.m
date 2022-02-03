function [NewPop] = Mallows_Ulam_sampling(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = Mallows_Ulam_sampling(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% Mallows_Ulam_sampling:  Samples a population of individuals from a Mallows Ulam model
%             
% INPUTS
% NumbVar:   Number of variables
% model:    model{1}  contains a vector of cumulated probabilities corresponding to the
%           the cumulated probabilities of having permutations at distance d
%             model{2} contains the consensus ranking.
%             model{3} contains calculated theta value.
% Card:      Vector with the dimension of all the variables. 
% AuxPop:    Auxiliary (selected) population (May be use for partial sampling or resampling)
% AuxFunVal: Evaluation of the data set (required for some sampling algorithms, not for this one)
% sampling_params{1}(1) = N: Number of generated individuals 
% OUTPUTS
% NewPop: Sampled individuals
%
% References:
% [1] C. L. Mallows: Non-null ranking models. Biometrika, 1957
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: Extending Distance-based Ranking Models in Estimation of Distribution Algorithms. In Proceedings of the Congress on Evolutionary Computation, 2014
%
% Created version 02/12/2015. Roberto Santana (roberto.santana@ehu.es) 
%
% Last version 04/22/2015. Roberto Santana (roberto.santana@ehu.es)  

m_probabilities = model{1}; % Distance cumulated probabilities matrix
m_consensus_ranking = model{2};
Theta = model{3};
m_problem_size = NumbVar;

N = cell2num(sampling_params{1}(1)); % Number of variables to sample

% 0.- Check if FerreShapes are stored 
  try
    eval(['load FerrerShapes_ncounts_',num2str(NumbVar),';']);   % The number of solutions at distance i will be in variable ncount
  catch
    disp(['File FerrerShapes_ncounts_',num2str(NumbVar),' is required to use this function. You can created using ComputeFerrerShapes.m ']); 
  end
  
  try
     eval(['load FerrerShapes_Lengths_',num2str(NumbVar),';']);   % The lengths of the shapes at each distance
  catch
    disp(['File FerrerShapes_Lengths_',num2str(NumbVar),' is required to use this function. You can created using ComputeFerrerShapes.m ']); 
  end

m_num_perms_distance = ncount(1:m_problem_size-1);                                   % Number of permutations at each distance [ncount is read from file FerrerShapes_ncounts]
%[N,cumsum(m_probabilities)]
ShapeIndices = sus(N,cumsum(m_probabilities));                                      % Determine the number of solutions to  sample at each distance  
for i=1:m_problem_size-1,                                                           % proportionally to the probabilities
  m_target_distances(i,:) = sum(ShapeIndices==i);
  if(m_target_distances(i,:)>0)
    aux_prob = FerrerShapes_Lengths{i}/sum(FerrerShapes_Lengths{i});
    AuxIndices = sus(m_target_distances(i,:),cumsum(aux_prob));
    split_samples_by_target_distance{i} = AuxIndices  ;    
  end 
end

 s=0;
 
 for target_distance=1:m_problem_size-1,
     
    if (m_target_distances(target_distance)~=0)          
          eval(['load permus_per_shape_',num2str(NumbVar),'_',num2str(target_distance),';']);         
          try
            eval(['load permus_per_shape_',num2str(NumbVar),'_',num2str(target_distance),';']);   % The lengths of the shapes at each distance
          catch
            disp(['File permus_per_shape_',num2str(NumbVar),'_',num2str(target_distance),' is required to use this function. You can created using ComputeFerrerShapes.m ']); 
          end        
          
       
          for j=1:length(FerrerShapes_Lengths{target_distance}),          
               shape = Ferrer_Shapes_l(j,:);                                    %  Select the corresponding Ferrer Shape
               In_Group_j = sum(split_samples_by_target_distance{target_distance}==j);        %  Determine how many solutions are sampled at this distance             
               for k=1:In_Group_j,                
                 %[j,k,s,shape  ]
                 m_composed = GeneratePermuAtLIS(NumbVar,shape);                %  Create a random permutation using Young Tables of compatible with Ferrer Shapes
                 m_sampling_permutation = Compose(m_composed,m_consensus_ranking);  % Create the permutation using the consensus ranking
                 s = s + 1;
                 NewPop(s,:)  = m_sampling_permutation;
                 %[target_distance,j,k,s]
               end
          end
    end     
 end
 
return;
