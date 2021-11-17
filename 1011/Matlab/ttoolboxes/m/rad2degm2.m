function val = rad2degm2(val)
% [ deg ] = rad2degm2( rad )
% Convert angles from radians to degrees, but does not change last row

    val(1:end-1,:) = (180/pi) * val(1:end-1,:);

end