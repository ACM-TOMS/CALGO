d = 0;
%For different degrees of Polya relaxation, change the value of d above.
for degP = 0:5
    clear Ai Bi Ci Di Pi dotPi mu LMIs T;

    A{1} = {[0 0 0],[0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0]};
    A{2} = {[1 0 0],[0 0 0 0; 0 0 0 0; -2 1 0 0; 0 0 0 0]};
    A{3} = {[1 0 1],[0 0 0 0; 0 0 0 0; 0 0 -1 0; 0 0 0 0]};
    A{4} = {[0 1 0],[0 0 0 0; 0 0 0 0; 0 0 0 0; 1 -1 0 0]};
    A{5} = {[0 1 1],[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 -1]};
    Ai = rolmipvar(A,'A',[2/3 2; 0.8 4/3; 1 3]);
    B{1} = {[1 0 0],[0;0;1;0]};
    Bi = rolmipvar(B,'B',[2/3 2]);
    Ci = rolmipvar([0 1 0 0],'C',0,0);
    Di = rolmipvar(0,'D',0,0);

    %%%%% Time-Invariant Parameters
    Pi = rolmipvar(4,4,'P','sym',[2 2 2],[degP degP degP]);
    mu = sdpvar(1,1);
    LMIs = [Pi > 0];
    T11 = Ai'*Pi + Pi*Ai + Ci'*Ci;
    T21 = Bi'*Pi + Di'*Ci;
    T22 = Di'*Di - mu*eye(1);
    T = [T11 T21'; T21 T22];
    %%%%%%%%
    
    T = polya(T,d);
    
    LMIs = [LMIs, T < 0];
    optimize(LMIs,mu);
    hinf(degP+1,d+1) = sqrt(double(mu))
end