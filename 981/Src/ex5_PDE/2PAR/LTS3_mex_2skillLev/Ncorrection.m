function NOPTScorr = Ncorrection(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,Ntol)

    %% correction to NOPTS for modified Talbot's method
    H = (Tmax-Tmin)/2;
    r = 18;
    N1max = NOPTS + H*(CONLAM*(CONNU+1)/2 + CONSIG/r);        % (4.4b)
    % initialization
    c  = (CONNU-1)/2;
    v  = (1 - (CONSIG-sigma0)/CONLAM)/2; % positive or negative
    r  = (1-v)/(1/2+c);
    y  = r*( 1/(exp(r)-1)-c );
    yp = (exp(r)*(1-r)-1)/(exp(r)-1)^2 - c;
    % iterations
    icount=0;
    while abs(y-v) > Ntol   &&   icount < 100
        icount = icount+1;
        r  = r + (v-y)/yp;
        y  = r*( 1/(exp(r)-1)-c );
        yp = (exp(r)*(1-r)-1)/(exp(r)-1)^2 - c;
    end
    N2max = NOPTS + H*((CONSIG+CONLAM)/r-CONLAM*(CONNU-1)/2); % (4.6b)
    NOPTScorr = round(max([N1max N2max]));
end

