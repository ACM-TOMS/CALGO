function checkRC(n, rc)

if length(n)==1
    if length(rc)==4
        if ~(rc(1)<=rc(2) && rc(3)<=rc(4) && rc(1)>0 ...
                && rc(3)>0 && rc(2)<=n && rc(4)<=n )
            error('Wrong rows/columns!')
        end
        
    elseif length(rc)==2
        if ~(rc(1)<=rc(2) && rc(1)>0  && rc(2)<=n)
            error('Wrong rows/columns!')
        end
        
    else
        error('Wrong Size!')
    end
else
    numBlk = length(n)-1;
    if ~(rc>0 && rc<=numBlk)
         error('Wrong block!')
    end
end
end