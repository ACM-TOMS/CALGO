function newPoly = affTransPoly(poly, A, b)
%%AFFTRANSPOLY   
%   This function applies the affine transformation 
%          x = A * y + b;
%   to the input polynomial f(x).
%
%   Usage:
%      newPoly = affTransPoly(poly, A, b);

    if ~checkPoly(poly)
        error('Input format is not correct');
    end

    At = A';
    calA = poly.supports';
    c = poly.coef;
    newPoly = [];

    % poly = \sum_{\alpha in calA} c_{\alpha}x^{\alpha}
    %      = \sum_{\alpha in calA} c_{\alpha} \Pi_{i=1}^n (At(:, i)' * y + b(i))^{\alpha_i}
    for kk = 1:size(calA, 2)
        termPoly.supports = sparse(1, size(At, 1));
        termPoly.coef = 1;
        % termPoly = c_{\alpha} x^{\alpha} 
        %          = c_{\alpha} \Pi_{i=1}^n (At(:, i)y + b(i))^{\alpha_i}

        alph = calA(:, kk); % \alpha
        nnzidx_alph = find(alph);
        for ii = nnzidx_alph(:)'
            varPoly = []; %(At(:, i)y + b(i))^{\alpha_i}

            Ati = At(:, ii);
            nnzidx_Ati = find(Ati);
            varPoly.supports = ...
                sparse(1:length(nnzidx_Ati), nnzidx_Ati, 1, length(nnzidx_Ati), length(Ati));
            varPoly.coef = Ati(nnzidx_Ati);
            if b(ii) ~= 0
                varPoly.supports = [varPoly.supports; sparse(1, length(Ati))];
                varPoly.coef = [varPoly.coef; b(ii)];
            end
            if alph(ii) > 1
                varPoly = powPoly(varPoly, alph(ii));
            end

            termPoly = multiplyPoly(varPoly, termPoly, false);
        end
        if ~isempty(termPoly)
            termPoly.coef = termPoly.coef * c(kk);
            newPoly = addPoly(newPoly, termPoly, false);
        end
    end
    %newPoly = simplifyPoly(newPoly);
end

