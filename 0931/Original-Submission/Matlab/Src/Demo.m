if verLessThan('matlab','7.9')
    warning('MATLAB:Ver',['MULTIPLICITY is designed for is designed for',...
        ' version 7.9 or higher (use in older version is your own risk).'])
end

disp('This is a master script to test all systems in the paper. ')
disp(['Please read the "computational experiments" section in the paper',...
    ' for detail.'])

addpath(pwd);
cd TestSuite

flag=waitinput(['Press enter to test cmbs1, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test cmbs1...')
Test_cmbs1();

flag=waitinput(['Press enter to test cmbs2, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test cmbs2...')
Test_cmbs2();

flag=waitinput(['Press enter to test mth191, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test mth191...')
Test_mth191();

flag=waitinput(['Press enter to test decker2, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test decker2...')
Test_decker2();

flag=waitinput(['Press enter to test Ojika2, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test Ojika2...')
Test_Ojika2(1);
Test_Ojika2(2);

flag=waitinput(['Press enter to test Ojika3, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test Ojika3...')
Test_Ojika3(1);
Test_Ojika3(2);

flag=waitinput(['Press enter to test Caprasse, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test Caprasse...')
Test_Caprasse();

flag=waitinput(['Press enter to test DZ1, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test DZ1...')
Test_DZ1();

flag=waitinput(['Press enter to test DZ2, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test DZ2...')
Test_DZ2();

flag=waitinput(['Press enter to test KSS, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test KSS...')
Test_KSS();

flag=waitinput(['Press enter to test cyclic cubic, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test cyclic cubic...')
Test_cyclic_cubic();

flag=waitinput(['Press enter to test ten-fold triangle, 0 to quit : ',...
    '(it will start automatically after 5 seconds) '],5);

if  ~isnan(flag) && ~flag
    return
end
disp(' ')
disp('Test ten-fold triangle...')
Test_ten_fold();

cd ..
