%from meibster at stackoverflow
function answer = isint(n)

if size(n) == [1 1]
    answer = isreal(n) && isnumeric(n) && round(n) == n &&  n >0;
else
    answer = false;
end


end