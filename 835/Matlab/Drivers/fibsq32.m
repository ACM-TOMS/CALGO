function [p,z] = fibsq32
%
%  test polynomial suggested by Goedecker
%  square of Fibocacci polynomial
%
    n = 32;
    p = [-1,ones(1,n)];
    p = conv(p,p);
    z = [-.9665892387165332, ... 
        -.9485092065867169-.1867222182129376*i, ... 
        -.9485092065867169+.1867222182129376*i, ... 
        -.8949223987510440-.3665688458356083*i, ... 
        -.8949223987510440+.3665688458356083*i, ... 
        -.8077651569138829-.5329094943253877*i, ... 
        -.8077651569138829+.5329094943253877*i, ... 
        -.6901871862102008-.6795950061744879*i, ... 
        -.6901871862102008+.6795950061744879*i, ... 
        -.5464383280299697-.8011754158225671*i, ... 
        -.5464383280299697+.8011754158225671*i, ... 
        -.3817161237254617-.8930914814237271*i, ... 
        -.3817161237254617+.8930914814237271*i, ... 
        -.2019801534661303-.9518321713734013*i, ... 
        -.2019801534661303+.9518321713734013*i, ... 
        -.1374081948052525e-1-.9750513807247308*i, ... 
        -.1374081948052525e-1+.9750513807247308*i, ... 
        .1761678077924798-.9616381291177512*i, ... 
        .1761678077924798+.9616381291177512*i, ... 
        .3608186423377573-.9117356125589683*i, ... 
        .3608186423377573+.9117356125589683*i, ... 
        .5334081767292508-.8267064135873833*i, ... 
        .5334081767292508+.8267064135873833*i, ... 
        .6874310990889808-.7090471025782057*i, ... 
        .6874310990889808+.7090471025782057*i, ... 
        .8167660514846921-.5622778345619681*i, ... 
        .8167660514846921+.5622778345619681*i, ... 
        .9156129648696007-.3909102172601782*i, ... 
        .9156129648696007+.3909102172601782*i, ... 
        .9783492503358520-.2007808902007197*i, ... 
        .9783492503358520+.2007808902007197*i, ... 
        1.999999999767169];
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
