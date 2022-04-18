%% runtests
% Environment to run tests of |IPscatt|.
%
% Warning: This functions uses |setInput|, such that a directory and
% several files are created. Also variables are cleared.
%
%% Syntax
%
%   seti = runtests()
%   seti = runtests(test)
%
%% Description
%
% * |seti = runtests()| runs all tests of |IPscatt| (same as |seti = runtests('all')|),
% i.e. tests derivative, mimo, compares near and far field results, tests grid
% scaling and adjoints.
% * |seti = runtests(test)| runs a specific test.
%
% *Specific tests by setting parameter |test|*
%
% * 'd'  : test derivative, see <testDerivative.html>.
% * 'm'  : test mimo (multiple input and multiple output), see <testMimo.html>.
% * 'mnf': test mimo (2D) and compare near and far field results, see <testMimoNearAndFar2D.html>.
% * 'g'  : test grid scaling, see <testGridScale.html>.
% * 'a'  : test adjoint of Kd and Kg, see <testAdjointKdAndKg.html>.
%
% * 'all' : runs all tests.
%
% *Figures*
%
% * Note that |inseti = 'tests'| is used, i.e. the parameters in 
% inseti/tests.m are used to influence the figure properties.
% * A folder |"date"T"time"_tests| is created, e.g.
% |20161018T154121_tests|, in folder |output| to save the created figures.
%
%% Examples
%
% *Example 1: Run all tests*
%
%   seti = runtests();
%
% *Example 2: Run the test of grid scaling*
%
%   seti = runtests('g');
%
%% Input Arguments
%
% * test    : string to choose a specific test, 
%             see Section Description for details
%             (set 'all' or empty to run all tests)
%
%% Output Arguments
%
% * seti    : structural array
%
%
% *Figures*
%
% (Fig. 1-10 are reserved because exp. set-up and predefined contrast is
% plotted, see <start.html>).
%
% Note that FF(q) does mean $\mathcal{F}(q)$.
%
% <html>
% <table>
% <tr><td><em>Fig. no.</em></td><td><em>Content</em></td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td><td><strong>01-10 Preliminary definitions...</strong><br> as described in <code>start.m</code> too</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td></td><td><code>expSetup.m</code></td></tr>
% <tr><td>01</td><td>source and measurement points plot...</td></tr>
% <tr><td></td><td><code>setContrast.m</code></td></tr>
% <tr><td>02</td><td>predefined contrast (real part)</td></tr>
% <tr><td>03</td><td>predefined contrast (imag part)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>04-10</td><td>currently unused</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td><td><strong>11-20 Test derivative (<code>testDerivative.m</code>)</strong><br>(2D: 11-14, 3D: 15-20)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>11</td><td>2D: test derivative of FF(q)</td></tr>
% <tr><td>12-14</td><td>currently unused</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>15</td><td>3D: test derivative of FF(q)</td></tr>
% <tr><td>16</td><td>3D (slice): Adjoint of derivative (ADFFq) computed by mimo</td></tr>
% <tr><td>17</td><td>3D (slice): Adjoint of derivative (ADFFq) computed by direct A,B decompostion</td></tr>
% <tr><td>18-20</td><td>currently unused</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td><td><strong>21-30 Test mimo (<code>testMimo.m</code>)</strong><br>(2D: 21-24, 3D: 25-30)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>21</td><td>2D: Computed scattered field</td></tr>
% <tr><td>22</td><td>2D: Reference scattered field</td></tr>
% <tr><td>23</td><td>2D: Difference</td></tr>
% <tr><td>24</td><td>currently unused</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>25</td><td>3D: Computed scattered field</td></tr>
% <tr><td>26</td><td>3D: Reference scattered field</td></tr>
% <tr><td>27</td><td>3D: Difference</td></tr>
% <tr><td>28</td><td>In 3D (slice): Computed scattered field</td></tr>
% <tr><td>29</td><td>In 3D (slice): Reference scattered field</td></tr>
% <tr><td>30</td><td>In 3D (slice): Difference</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td><td><strong>31-40 Test grid scaling (<code>testGridScale.m</code>)</strong><br>In 2D case: (in 3D no plot output)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>31</td><td>original matrix A (fine grid)</td></tr>
% <tr><td>32</td><td>matrix A down scaled</td></tr>
% <tr><td>33</td><td>matrix A after down and up scale</td></tr>
% <tr><td>34-40</td><td>currently unused</td></tr>
% </table>
% </html>
%
%% More About
%
% * In diary.log in the created folder inside the folder output you can
% find the terminal output.
% * Not all results are plotted (so look in terminal or in diary.log).
% * Not all plots are saved (e.g. in <testMimo.html> for different nCD, only the
% last one is saved)
%
% *How to run tests standalone?*
%
% Call before:
%
%   init;
%   setInput;
%
% Better: Use runtests, e.g.:
%
%   seti = runtests('g');
%
%% See Also
%
% * <tests.html>
% * <testDerivative.html>
% * <testMimo.html>
% * <testMimoNearAndFar2D.html>
% * <testGridScale.html>
% * <testAdjointKdAndKg.html>
%
%% Code
%
function seti = runtests(test)

disp(' ')
disp('################ -------- tests of IPscatt ------- ################')
disp(' ')

inseti = 'tests'; % the variable inseti will be used in the routine setInput

%%
init;
setInput; % Clears all variables (except closed and seti).
          % This is not recommended in tests... but for standalone running it is necessary...

%% Diary log: start
disp(' ')
disp('## diary -- start diary log ---------------------------------------')
disp(' ')
diary(sprintf('%s/diary%s.log',seti.dirname,seti.fileSuffix)); % store all output in diary
% The structure array seti was defined before in setInput.
diary on; % later diary off

%% Run options

enablePauses = 0;

%%

if nargin == 0
    test = 'all';
end

switch test
        %% Tests
    case 'd'
        disp(' ')
        disp('## testDerivative -- test derivative ----------------------')
        disp(' ')
        testDerivative;
    case 'm'
        disp(' ')
        disp('## testMimo -- test multiple input and multiple output ----')
        disp(' ')
        testMimo;
    case 'mnf'
        disp(' ')
        disp('## testMimoNearAndFar2D -- test mimo (2D) and compare near and far field results --')
        disp(' ')
        testMimoNearAndFar2D;
    case 'g'
        disp(' ')
        disp('## testGridScale -- testGridScale -------------------------')
        disp(' ')
        testGridScale;
    case 'a'
        disp(' ')
        disp('## testAdjointKdAndKg -- test adjoint of Kd and Kg --------')
        disp(' ')
        testAdjointKdAndKg;
    case 'all' % all tests

        disp(' ')
        disp('## tests: all -- run all tests-----------------------------')
        disp(' ')
        %% All tests: tests with fig output
        
        disp(' ')
        disp('## testDerivative -- 2D -----------------------------------')
        disp(' ')
        seti.dim = 2;
        testDerivative;
        disp(' '); disp('-----'); disp(' ')
        drawnow;
        if enablePauses
            disp('(Press any key to continue)')
            pause
        end
        
        disp(' ')
        disp('## testDerivative -- 3D -----------------------------------')
        disp(' ')
        seti.dim = 3;
        seti.incPntsType = 'sphereFibo';
        seti.measPntsType = 'sphereFibo';
        testDerivative;
        seti = rmfield(seti,'incPntsType');
        seti = rmfield(seti,'measPntsType');
        disp(' '); disp('-----'); disp(' ')
        drawnow;
        if enablePauses
            disp('(Press any key to continue)')
            pause
        end
        
        disp(' ')
        disp('## testMimo -- 2D -----------------------------------------')
        disp(' ')
        seti.dim = 2;
        seti.model = 'helmholtz2D';
        testMimo;
        disp(' '); disp('-----'); disp(' ')
        drawnow;
        if enablePauses
            disp('(Press any key to continue)')
            pause
        end

        disp(' ')
        disp('## testMimo -- 3D -----------------------------------------')
        disp(' ')
        seti.dim = 3;
        seti.model = 'helmholtz3D';
        seti.incPntsType = 'sphereFibo';
        seti.measPntsType = 'sphereFibo';
        testMimo;
        seti = rmfield(seti,'incPntsType');
        seti = rmfield(seti,'measPntsType');
        disp(' '); disp('-----'); disp(' ')
        drawnow;
        if enablePauses
            disp('(Press any key to continue)')
            pause
        end

        disp(' ')
        disp('## testGridscale -- 2D ------------------------------------')
        disp(' ')
        dim = 2; % The parameter dim is used in testGridScale.
        testGridScale;
        disp(' '); disp('-----'); disp(' ')
        drawnow;
        if enablePauses
            disp('(Press any key to continue)')
            pause
        end
        
        disp(' ')
        disp('## testGridScale -- 3D ------------------------------------')
        disp(' ')
        dim = 3;
        testGridScale; % The parameter dim is used in testGridScale.
        disp(' '); disp('-----'); disp(' ')
        drawnow;
        if enablePauses
            disp('(Press any key to continue)')
            pause
        end

        %% All tests: tests without fig output

        disp(' ')
        disp('## Test mimo (2D) and compare near and far field results --')
        disp(' ')
        testMimoNearAndFar2D;
        disp(' '); disp('-----'); disp(' ')

        disp(' ')
        disp('## Test adjoint of Kd and Kg ------------------------------')
        disp(' ')
        testAdjointKdAndKg;
        disp(' '); disp('-----'); disp(' ')

end

%% Diary log: end
disp(' ')
disp('## diary -- end diary log -----------------------------------------')
disp(' ')
diary off;

end
