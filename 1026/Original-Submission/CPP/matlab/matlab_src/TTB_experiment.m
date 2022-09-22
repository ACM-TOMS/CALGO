% Path to TTB
home = getenv('HOME');
addpath(home + "/projects/tensor_toolbox/");
addpath(home + "/projects/CP-CALS-priv/data/");

n = maxNumCompThreads();
if n ~= 1
  n = maxNumCompThreads(str2num(getenv('OMP_NUM_THREADS')));
end
n = maxNumCompThreads();
%if n == 24
%	maxNumCompThreads(12)
%	n = maxNumCompThreads
%	if n ~= 12
%		disp("error with number of threads")
%		exit
%	end
if n == 1
	disp("Threads = 1")
elseif n == 24
  disp("Threads = 24")
else
  fprintf('Unrecognized number of threads %d\n', n);
	exit
end

clc
clear all
echo off
rand('state', 0)

rank_min = 1;
rank_max = 20;
copies = 20;

%% EEMData
tic;
load('eemdata.mat')
X = fluor.data(:,:,:,1);
sizes = size(X);
mode1 = sizes(1);
mode2 = sizes(2);
mode3 = sizes(3);

X(isnan(X))=0;
X = tensor(X);
% X = tensor(rand(mode1, mode2, mode3));
%X = full(ktensor({rand(mode1, tensor_rank), rand(mode2, tensor_rank), rand(mode3, tensor_rank)}));

ranks = zeros(1, 400);
k = 1;
for i = rank_min:rank_max
    for j = 1:copies
        ranks(k) = i;
        k = k + 1;
    end
end
rank_sum = sum(ranks);
fprintf('Ktensors: %d\n', size(ranks, 2));
fprintf('Ktensor rank sum: %d\n', rank_sum);

ktensors = {};
for r = 1:size(ranks, 2)
    ktensors{r} = ktensor({rand(mode1, ranks(r)),
                           rand(mode2, ranks(r)),
                           rand(mode3, ranks(r))});
end
time_gen = toc;
fprintf('Time to generate tensor and models: %f\n', time_gen);

M = {};
disp(' ')
disp('START TTB ALS')

timings = zeros(1, size(ktensors, 2));
for k = 1:size(ktensors, 2)
	tic
    M{k} = cp_als(X, size(ktensors{k}.lambda, 1), 'maxiters', 1000, 'tol', 1e-6, 'init', ktensors{k}.U);
    disp(k);
 	timings(k)= toc;
end
path = "../../data/MKL/TTB_" + string(mode1) + "-" + string(mode2) + "-" + string(mode3) + "_" + string(maxNumCompThreads) + ".csv"
fileID = fopen(path, 'w');
fprintf(fileID, '%10f\n', timings);

disp('STOP TTB ALS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
echo off
rand('state', 0)

rank_min = 1;
rank_max = 20;
copies = 20;

%% Cancer Data
tic;
load('fluorescence_data.mat')
X = X_UD.data;
sizes = size(X);
mode1 = sizes(1);
mode2 = sizes(2);
mode3 = sizes(3);

X(isnan(X))=0;
X = tensor(X);
% X = tensor(rand(mode1, mode2, mode3));
%X = full(ktensor({rand(mode1, tensor_rank), rand(mode2, tensor_rank), rand(mode3, tensor_rank)}));

ranks = zeros(1, 400);
k = 1;
for i = rank_min:rank_max
    for j = 1:copies
        ranks(k) = i;
        k = k + 1;
    end
end
rank_sum = sum(ranks);
fprintf('Ktensors: %d\n', size(ranks, 2));
fprintf('Ktensor rank sum: %d\n', rank_sum);

ktensors = {};
for r = 1:size(ranks, 2)
    ktensors{r} = ktensor({rand(mode1, ranks(r)),
                           rand(mode2, ranks(r)),
                           rand(mode3, ranks(r))});
end
time_gen = toc;
fprintf('Time to generate tensor and models: %f\n', time_gen);

M = {};
disp(' ')
disp('START TTB ALS')

timings = zeros(1, size(ktensors, 2));
for k = 1:size(ktensors, 2)
	tic
    M{k} = cp_als(X, size(ktensors{k}.lambda, 1), 'maxiters', 1000, 'tol', 1e-6, 'init', ktensors{k}.U);
 	timings(k)= toc;
end
path = "../../data/MKL/TTB_" + string(mode1) + "-" + string(mode2) + "-" + string(mode3) + "_" + string(maxNumCompThreads) + ".csv"
fileID = fopen(path, 'w');
fprintf(fileID, '%10f\n', timings);

disp('STOP TTB ALS')
