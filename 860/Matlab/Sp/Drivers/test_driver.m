% Create matrix with the 2-D simplex vertices as colunms
S2 = [0  1.0  0.6
     0  0.2  1.0]
% Compute the 10^2 subdivision
T=simples(10,2);
% Draw the subdivision
draw_simplices(S2,T)

% Create matrix with the 3-D simplex vertices as colunms
S3 = [0.0    0.5    1.0    0.5
     1.0    0.0    1.0    0.5
     0.0    0.0    0.0    1.0]
% Compute the 6^3 subdivision
T=simples(6,3);
% Draw the subdivision
draw_simplices(S3,T)
view(-15,15)


