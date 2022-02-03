function s = order2string(order)

s = ' ';
if (order<=4)
    for i=1:order
        s = strcat(s,'''');
    end
else
    s = strcat('^(', num2str(order), ')');
end
end