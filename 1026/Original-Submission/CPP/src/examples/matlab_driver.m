% Path to your CALS MEX
addpath("/home/chris/projects/SimALS/cmake-build-matlab/")
addpath("/home/chris/projects/SimALS/matlab/matlab_src/")
% Path to TTB
addpath("/home/chris/projects/tensor_toolbox/")

clc
clear all
echo off
rand('state', 0)

mode1 = 210;
mode2 = 210;
mode3 = 210;
rank_min = 1;
rank_max = 5;
copies = 10;

tic;
X = full(ktensor({rand(mode1, 5), rand(mode2, 5), rand(mode3, 5)}));
ranks = zeros(1, 50);
k = 1;
for i = rank_min:rank_max
    for j = 1:copies
        ranks(k) = i;
        k = k + 1;
    end    
end
rank_sum = sum(ranks);

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
disp('TTB ALS...')
tic
for k = 1:size(ktensors, 2)
    M{k} = cp_als(X, size(ktensors{k}.lambda, 1), 'maxiters', 50, 'tol', 1e-4, 'init', ktensors{k}.U, 'printitn', 0);
end
time_tt = toc;
pause(3)
disp(' ')
disp('CALS...')
tic
M1 = cp_cals(X, ranks, ktensors, 'mttkrp-method', 'auto', 'tol', 1e-4, 'maxiters', 50, 'buffer-size',  rank_sum, 'no-cuda');
time_cals = toc;

t = size(M, 2) * [];
for i = 1:size(M,2)
    t(i) = abs(norm(X) - norm(M{i})) - abs(norm(X) - norm(M1{i}));
end
disp(' ')
disp('-----------------------------------------------------------')
fprintf('ALS(TTB) time: %0.5f\n', time_tt);
fprintf('CALS time    : %0.5f\n', time_cals);
fprintf('Speedup      : %0.2f\n', time_tt / time_cals);

