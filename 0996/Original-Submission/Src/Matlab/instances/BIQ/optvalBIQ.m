function optval = optvalBIQ(instance)

    if nargin == 0
        optval = [];
    end

    k = strfind(instance,'-');
    if k > 0
        instance = [instance(1:k-1),'_',instance(k+1:end)];
    end

    optval.bqp100_1 = -7970;
    optval.bqp100_2 = -11036;
    optval.bqp100_3 = -12723;
    optval.bqp100_4 = -10368;
    optval.bqp100_5 = -9083;
    optval.bqp100_6 = -10210;
    optval.bqp100_7 = -10125;
    optval.bqp100_8 = -11435;
    optval.bqp100_9 = -11455;
    optval.bqp100_10 = -12565;
    optval.bqp250_1 = -45607;
    optval.bqp250_2 = -44810;
    optval.bqp250_3 = -49037;
    optval.bqp250_4 = -41274;
    optval.bqp250_5 = -47961;
    optval.bqp250_6 = -41014;
    optval.bqp250_7 = -46757;
    optval.bqp250_8 = -35726;
    optval.bqp250_9 = -48916;
    optval.bqp250_10 = -40442;
    optval.bqp500_1 = -116586;
    optval.bqp500_2 = -128223;
    optval.bqp500_3 = -130812;
    optval.bqp500_4 = -130097;
    optval.bqp500_5 = -125487;
    optval.bqp500_6 = -121772;
    optval.bqp500_7 = -122201;
    optval.bqp500_8 = -123559;
    optval.bqp500_9 = -120798;
    optval.bqp500_10 = -130619;
    optval.gka1a = -3414;
    optval.gka2a = -6063;
    optval.gka3a = -6037;
    optval.gka4a = -8598;
    optval.gka5a = -5737;
    optval.gka6a = -3980;
    optval.gka7a = -4541;
    optval.gka8a = -11109;
    optval.gka1b = -133;
    optval.gka2b = -121;
    optval.gka3b = -118;
    optval.gka4b = -129;
    optval.gka5b = -150;
    optval.gka6b = -146;
    optval.gka7b = -160;
    optval.gka8b = -145;
    optval.gka9b = -137;
    optval.gka10b = -154;
    optval.gka1c = -5058;
    optval.gka2c = -6213;
    optval.gka3c = -6665;
    optval.gka4c = -7398;
    optval.gka5c = -7362;
    optval.gka6c = -5824;
    optval.gka7c = -7225;
    optval.gka1d = -6333;
    optval.gka2d = -6579;
    optval.gka3d = -9261;
    optval.gka4d = -10727;
    optval.gka5d = -11626;
    optval.gka6d = -14207;
    optval.gka7d = -14476;
    optval.gka8d = -16352;
    optval.gka9d = -15656;
    optval.gka10d = -19102;
    optval.gka1e = -16464;
    optval.gka2e = -23395;
    optval.gka3e = -25243;
    optval.gka4e = -35594;
    optval.gka5e = -35154;
    optval.gka1f = -61194;
    optval.gka2f = -100161;
    optval.gka3f = -138035;
    optval.gka4f = -172771;
    optval.gka5f = -190507;

    if isfield(optval, instance)
	    optval = optval.(instance);
    else
	    optval = 1e10;
    end

end

