function [objPoly, I01, Icomp, relaxOrder, params] = example1
%% EXAMPLE1
%
% Usage: [objPoly, I01, Icomp, relaxOrder, params] = EXAMPLE1
%
% This function outputs the sparsePOP representation of the following POP:
%
% % \begin{equation}
% % \begin{array}{ll}
% %    \mbox{minimize} & f_0(\x) = 0.5x_1 -1.8x_3 -2.2x_5 +3x_3^2 +x_1x_2x_4 + 1.3x_2x_4x_5      \\
% %    \mbox{subject to} &  x_2x_3 = 0,\  x_3x_4 = 0, \  x_1, x_2 \in \{0,1\}, \  x_3, x_4, x_5 \in [0,1]  
% % \end{array}
% % \end{equation}
%
    objPoly.typeCone   = 1;
    objPoly.sizecone   = 1;
    objPoly.dimVar     = 5;
    objPoly.degree     = 3;
    objPoly.noTerms    = 6;
    objPoly.supports   = ...
        [1 0 0 0 0;
        0 0 1 0 0;
        0 0 0 0 1;
        0 0 2 0 0;
        1 1 0 1 0;
        0 1 0 1 1];
    objPoly.coef = [0.5; -1.8; -2.2; 3; 1; 1.3];
    Icomp = logical([0 1 1 0 0; 0 0 1 1 0]);
    I01 = logical([1 1 0 0 0]);
    
    relaxOrder = 2;
    params.sparseSW = 1;
end
