function val = deg2radm2(val)
% [ rad ] = deg2radm2( deg )
% Convert angles from degrees to radians. Does not convert the last row.

val(1:end-1,:) = (pi/180) * val(1:end-1,:);

end