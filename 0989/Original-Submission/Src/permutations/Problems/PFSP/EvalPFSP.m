function[val] =  EvalPFSP(vector)
% [val] =  EvalPFSP(vector)
% EvalPFSP: Using a processing time matrix computes the cost of executing a given permutation of jobs.
% Input:
% vector: Permutation vector
% val:  Cost of the permutation according to global PFSPInstance processingtimes variable
%
% NOTE: In order to use the function correctly you must declare a global
% boolean variable depending on these cases:
%       - true: The permutation will be inversed before the evaluation.
%       - false: The permutation will NOT be inversed.
%
%   By default the value is true
%
% Created version 12/19/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 03/19/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es)

global PFSPInstance;
global inverse;

%if inverse is not defined by the user, put the default value.
if isempty(inverse)
    inverse=true;
end

PermutatioSize = size(vector,2);
val = 0;

[processingtimes, machines, jobs] = PFSPInstance{:}; % Get the variables from the processed instance
timeTable = zeros(1,machines);

if inverse
	genes=Invert(vector);
else
	genes=vector;
end

prev_machine=0;
first_gene=genes(1);

j=1:machines;
timeTable = cumsum(processingtimes(j,first_gene));

%Calculate the fitness by simulating the execution of a given permutation of jobs
fitness = timeTable(machines);
for z=2:jobs
	job = genes(z);

	timeTable(1) =timeTable(1) + processingtimes(1,job);
	prev_machine=timeTable(1);

	for machine=2:machines
		timeTable(machine) = max(prev_machine,timeTable(machine))+ processingtimes(machine,job);
		prev_machine=timeTable(machine);
	end

	fitness = fitness + timeTable(machines);
end
	
val = -fitness;