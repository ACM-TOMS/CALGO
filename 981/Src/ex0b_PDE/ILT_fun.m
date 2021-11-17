function u = ILT_fun( eta,t )
    % return u(eta,t), the true solution

    % phi(eta,t), eta>0, t>0
    %beta=2; a=1;
    phi = @(n,x) 1/2*(exp(-n).*erfc(n./(2*sqrt(x))-sqrt(x)) + exp(+n).*erfc(n./(2*sqrt(x))+sqrt(x)));
    % f'(t)
    f1 = @(x) sin(3*x)/6 + x.*cos(3*x)/2;
    %integrand = @(eta,t,tau) g(t-tau).*Phi(eta,tau);
    
    %% u as convolution of g and Phi, in [0,t]
     u = exp(-eta).*quadgk( @(tau) f1(t-tau).*phi(eta,tau), 0,t,'AbsTol',1e-12 );

end
