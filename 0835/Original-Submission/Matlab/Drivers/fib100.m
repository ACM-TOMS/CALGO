function [p,z] = fib100
%
%  test polynomial suggested by Goedecker
%
    n = 100;
    p = [-1,ones(1,n)];
    z = [-.9891099736226092, ... 
        -.9871753104718534-.6190245872746059e-1*i, ... 
        -.9871753104718534+.6190245872746059e-1*i, ... 
        -.9813788056446001-.1235640063982463*i, ... 
        -.9813788056446001+.1235640063982463*i, ... 
        -.9717428840956687-.1847446602565987*i, ... 
        -.9717428840956687+.1847446602565987*i, ... 
        -.9583048244621342-.2452062905225305*i, ... 
        -.9583048244621342+.2452062905225305*i, ... 
        -.9411166150274026-.3047135378280153*i, ... 
        -.9411166150274026+.3047135378280153*i, ... 
        -.9202447528605484-.3630347198287647*i, ... 
        -.9202447528605484+.3630347198287647*i, ... 
        -.8957699869173070-.4199427234459649*i, ... 
        -.8957699869173070+.4199427234459649*i, ... 
        -.8677870061092108-.4752158792455410*i, ... 
        -.8677870061092108+.4752158792455410*i, ... 
        -.8364040735640849-.5286388145285358*i, ... 
        -.8364040735640849+.5286388145285358*i, ... 
        -.8017426085137990-.5800032817847029*i, ... 
        -.8017426085137990+.5800032817847029*i, ... 
        -.7639367174532022-.6291089592520593*i, ... 
        -.7639367174532022+.6291089592520593*i, ... 
        -.7231326764169818-.6757642204274601*i, ... 
        -.7231326764169818+.6757642204274601*i, ... 
        -.6794883664183336-.7197868694867672*i, ... 
        -.6794883664183336+.7197868694867672*i, ... 
        -.6331726642844369-.7610048396973263*i, ... 
        -.6331726642844369+.7610048396973263*i, ... 
        -.5843647913085710-.7992568520396451*i, ... 
        -.5843647913085710+.7992568520396451*i, ... 
        -.5332536223172150-.8343930313987525*i, ... 
        -.5332536223172150+.8343930313987525*i, ... 
        -.4800369579227939-.8662754778380540*i, ... 
        -.4800369579227939+.8662754778380540*i, ... 
        -.4249207628992477-.8947787906289414*i, ... 
        -.4249207628992477+.8947787906289414*i, ... 
        -.3681183737790589-.9197905428773492*i, ... 
        -.3681183737790589+.9197905428773492*i, ... 
        -.3098496789278928-.9412117047633791*i, ... 
        -.3098496789278928+.9412117047633791*i, ... 
        -.2503402745082667-.9589570135917125*i, ... 
        -.2503402745082667+.9589570135917125*i, ... 
        -.1898205998989927-.9729552890388070*i, ... 
        -.1898205998989927+.9729552890388070*i, ... 
        -.1285250562957008-.9831496921783526*i, ... 
        -.1285250562957008+.9831496921783526*i, ... 
        -.6669111238377778e-1-.9894979270705116*i, ... 
        -.6669111238377778e-1+.9894979270705116*i, ... 
        -.4558401154059704e-2-.9919723839157337*i, ... 
        -.4558401154059704e-2+.9919723839157337*i, ... 
        .5763218786929409e-1-.9905602230050737*i, ... 
        .5763218786929409e-1+.9905602230050737*i, ... 
        .1196394164912155-.9852633989536015*i, ... 
        .1196394164912155+.9852633989536015*i, ... 
        .1812226010328281-.9760986249939283*i, ... 
        .1812226010328281+.9760986249939283*i, ... 
        .2421425057566853-.9630972774521361*i, ... 
        .2421425057566853+.9630972774521361*i, ... 
        .3021622210555556-.9463052409577104*i, ... 
        .3021622210555556+.9463052409577104*i, ... 
        .3610480204611693-.9257826954965046*i, ... 
        .3610480204611693+.9257826954965046*i, ... 
        .4185701899698170-.9016038471668328*i, ... 
        .4185701899698170+.9016038471668328*i, ... 
        .4745038225339292-.8738566055392182*i, ... 
        .4745038225339292+.8738566055392182*i, ... 
        .5286295698615102-.8426422119868863*i, ... 
        .5286295698615102+.8426422119868863*i, ... 
        .5807343429795305-.8080748254361926*i, ... 
        .5807343429795305+.8080748254361926*i, ... 
        .6306119525307706-.7702810749337576*i, ... 
        .6306119525307706+.7702810749337576*i, ... 
        .6780636798327496-.7293995925445015*i, ... 
        .6780636798327496+.7293995925445015*i, ... 
        .7228987709578736-.6855805456977909*i, ... 
        .7228987709578736+.6855805456977909*i, ... 
        .7649348495515867-.6389851953915993*i, ... 
        .7649348495515867+.6389851953915993*i, ... 
        .8039982514449404-.5897855154613771*i, ... 
        .8039982514449404+.5897855154613771*i, ... 
        .8399242976881979-.5381639173058062*i, ... 
        .8399242976881979+.5381639173058062*i, ... 
        .8725575452576575-.4843131310830981*i, ... 
        .8725575452576575+.4843131310830981*i, ... 
        .9017520886609776-.4284362924194566*i, ... 
        .9017520886609776+.4284362924194566*i, ... 
        .9273720302113955-.3707472629394290*i, ... 
        .9273720302113955+.3707472629394290*i, ... 
        .9492922833419482-.3114711598054638*i, ... 
        .9492922833419482+.3114711598054638*i, ... 
        .9673998997347004-.2508449722447981*i, ... 
        .9673998997347004+.2508449722447981*i, ... 
        .9815960786437149-.1891180053506043*i, ... 
        .9815960786437149+.1891180053506043*i, ... 
        .9917988816037989-.1265517518780266*i, ... 
        .9917988816037989+.1265517518780266*i, ... 
        .9979464229745987-.6341873580445710e-1*i, ... 
        .9979464229745987+.6341873580445710e-1*i, ... 
        2.000000000000000];
    z = [z',ones(n,1)];
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
