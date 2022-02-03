function [closes] = item_closes_cycle(pos, item, permu)

while(item < pos)
    item = permu(item);
end
if(item == pos)
    closes=true;
    return
end
closes=false;
return



