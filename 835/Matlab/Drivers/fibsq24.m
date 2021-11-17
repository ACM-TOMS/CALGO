function [p,z] = fibsq24
%
%  test polynomial suggested by Goedecker
%  square of Fibocacci polynomial
%
    n = 24;
    p = [-1,ones(1,n)];
    p = conv(p,p);
    z = [-.9558467190743565, ... 
        -.9244083656460181-.2442575811931259*i, ... 
        -.9244083656460181+.2442575811931259*i, ... 
        -.8320677914606713-.4727784060949767*i, ... 
        -.8320677914606713+.4727784060949767*i, ... 
        -.6846250033032301-.6707962641470166*i, ... 
        -.6846250033032301+.6707962641470166*i, ... 
        -.4913433576989534-.8254218906657225*i, ... 
        -.4913433576989534+.8254218906657225*i, ... 
        -.2643730603203030-.9264240880183555*i, ... 
        -.2643730603203030+.9264240880183555*i, ... 
        -.1800084671534713e-1-.9668298465915332*i, ... 
        -.1800084671534713e-1+.9668298465915332*i, ... 
        .2322158170552186-.9432943625623267*i, ... 
        .2322158170552186+.9432943625623267*i, ... 
        .4703449187206489-.8561989307649443*i, ... 
        .4703449187206489+.8561989307649443*i, ... 
        .6808619003199168-.7094507040001145*i, ... 
        .6808619003199168+.7094507040001145*i, ... 
        .8490194865023282-.5100588510854889*i, ... 
        .8490194865023282+.5100588510854889*i, ... 
        .9602996918859324-.2681678908797330*i, ... 
        .9602996918859324+.2681678908797330*i, ... 
        1.999999940395313];
    z = [z',2*ones(n,1)];
    if norm(imag(z(:,1))) == 0 
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.15f \t \t \t %3g \n', z');
    else
        fprintf('                 roots ')
        fprintf('   \t\t\t\t\t\t     multiplicities\n');
        fprintf('\n');
        fprintf('%22.15f + %22.15f i \t \t %3g \n', ...
            [real(z(:,1)),imag(z(:,1)),z(:,2)]');

    end;    
