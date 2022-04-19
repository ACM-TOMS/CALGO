function [NewPop] = InitPermutations(NumbVar,PopSize,Card,seeding_params)
% [NewPop] = InitPermutations(NumbVar,PopSize,Card,sampling_params)
% InitPermutations:         Random Initialization of a population
%                           comprising permutations
% INPUTS
% NumbVar: Number of variables
% Card: For discrete variables:    Vector with the dimension of all the variables. 
% OUTPUTS
% NewPop: Sampled individuals
%


filename='poblacion PFSP prueba.txt';
%filename='poblacion LOP prueba.txt';
%filename='poblacion QAP prueba.txt';
%filename='poblacion TSP prueba.txt';
%Read the instance file
fid = fopen(filename);

%%Read first line
%tline = fgetl(fid);


%Read the variables
% fscanf fills the array in column order
[B, count] = fscanf(fid,'%d',[NumbVar,PopSize]); %general take [NumbVar,PopSize] pop %FPSP = [20,200], %LOP= [50,500]
fclose(fid);
NewPop = B' + 1; % +1 para usar indices de 1..N


global randMat;

%filename='numeros random seleccion PFSP.txt';
%filename='numeros random seleccion LOP.txt';
%filename='numeros random seleccion QAP.txt';
%filename='numeros random seleccion TSP.txt';

filename='numeros random seleccion C PFSP.txt';
%filename='numeros random seleccion C LOP.txt';
%filename='numeros random seleccion C QAP.txt';
%filename='numeros random seleccion C TSP.txt';

%filename='numeros random seleccion GC PFSP.txt';
%filename='numeros random seleccion GC LOP.txt';
%filename='numeros random seleccion GC QAP.txt';
%filename='numeros random seleccion GC TSP.txt';
%filename='numeros random seleccion  best perm GC TSP.txt';

%Read the instance file
fid = fopen(filename);
%Read the variables
% fscanf fills the array in column order
randMat = zeros(PopSize-1,NumbVar-1);
for i=1:PopSize-1
    %Read first line
    tline = fgetl(fid);
    [B, count] = fscanf(fid,'%f',[1,NumbVar-1]);
    randMat(i,:) = B';
    %Read end of line
    tline = fgetl(fid);
end
fclose(fid);

 
