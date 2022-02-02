function ineqPolySys = lin2Poly(A, b)
%%% LIN2POLY
%  Ax = b  ->  sparsePOP format
%  Usage:  eqPolySys = LIN2POLY(A, b);

    b = b(:);
    [mm, nn] = size(A);
    if length(b) ~= mm
        error('The size of A and b is inconsistent');
    elseif nn == 0
        error('No variable');
    end

    ineqPolySys = cell(mm, 1);
    
    At = A'; % for acceleration
    for ii = 1:mm
        %ineqPolySys{ii}.typeCone = -1;
        %ineqPolySys{ii}.supports = [sparse(nn, 1); speye(nn)];
        ineqPolySys{ii}.supports = [sparse(1, nn); speye(nn)];
        ineqPolySys{ii}.coef = [b(ii); At(:, ii)];
    end
end