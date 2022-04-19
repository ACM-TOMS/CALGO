classdef SolverTest < matlab.unittest.TestCase
    % SolverTest tests solutions to the functions of the multicomplex class
    % definition
    
    methods (Test)
        function testmatrep(testCase)
            actSolution = matrep(multicomplex([1,2,3,4]));
            expSolution = [1    -2    -3     4
                           2     1    -4    -3
                           3    -4     1    -2
                           4     3     2     1];
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testplus(testCase)
            actSolution = multicomplex([1,-2])+multicomplex([-2,-2]);
            expSolution = multicomplex([-1,-4]);
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testminus(testCase)
            actSolution = multicomplex([1,-2])-multicomplex([-2,-2]);
            expSolution = multicomplex([3,0]);
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testless_than(testCase)
            actSolution = multicomplex([1,-2])<multicomplex([-2,-2]);
            expSolution = false;
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testmore_than(testCase)
            actSolution = multicomplex([1,-2])>multicomplex([-2,-2]);
            expSolution = true;
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testless_than_eq_to(testCase)
            actSolution = multicomplex([1,-2])<=multicomplex([-2,-2]);
            expSolution = false;
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testmore_than_eq_to(testCase)
            actSolution = multicomplex([1,-2])>=multicomplex([-2,-2]);
            expSolution = true;
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testequal(testCase)
            actSolution = multicomplex([1,-2])==multicomplex([1,-2]);
            expSolution = true;
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testnot_equal(testCase)
            actSolution = multicomplex([1,-2])~=multicomplex([1,-2]);
            expSolution = false;
            testCase.verifyEqual(actSolution,expSolution);
        end
        
        
        function testmultiplication(testCase)
            actSolution = multicomplex([1,2,3,4])*multicomplex([-2,-4,-2.5,-1]);
            expSolution = multicomplex([9.5000 5 9.5000 -26]);
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testdivision(testCase)
            actSolution = imgn(multicomplex([1,2,3,4])/multicomplex([-2,-4,-2.5,-1]));
            expSolution = -0.258959537572254;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testpower(testCase)
            actSolution = imgn(multicomplex([-2,-4,-2.5,-1])^-0.45);
            expSolution = -0.167913799065957;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        
        
        function testatan2(testCase)
            actSolution = imgn(atan2(multicomplex([-2,-4,-2.5,-1]),multicomplex([3,-2])));
            expSolution = 0.574946838860950;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testatan(testCase)
            actSolution = imgn(atan(multicomplex([-2,-4,-2.5,-1])));
            expSolution = 0.120571907821568;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testlog(testCase)
            actSolution = imgn(log(multicomplex([-2,-4,-2.5,-1])));
            expSolution = -0.336657276181865;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testsin(testCase)
            actSolution = imgn(sin(multicomplex([-2,-4,-2.5,-1])));
            expSolution = 139.758676491826;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testcos(testCase)
            actSolution = imgn(cos(multicomplex([-2,-4,-2.5,-1])));
            expSolution = -91.0090265276186;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testsinh(testCase)
            actSolution = imgn(sinh(multicomplex([-2,-4,-2.5,-1])));
            expSolution = 0.219525252215660;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testcosh(testCase)
            actSolution = imgn(cosh(multicomplex([-2,-4,-2.5,-1])));
            expSolution = -0.397397608271755;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testtan(testCase)
            actSolution = imgn(tan(multicomplex([-2,-4,-2.5,-1])));
            expSolution = 0.012670705465629;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testasin(testCase)
            actSolution = imgn(asin(multicomplex([0.5,-4,-2.5,-1])));
            expSolution = -0.248930526780233;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testacos(testCase)
            actSolution = imgn(acos(multicomplex([0.5,-4,-2.5,-1])));
            expSolution = 0.248930526780237;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testexp(testCase)
            actSolution = imgn(exp(multicomplex([-2,-4,-2.5,-1])));
            expSolution = -0.177872356056095;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        
        
        
        function testlast_img_coeff(testCase)
            actSolution = imgn(multicomplex([1,2,3,4]));
            expSolution = 4;
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testCX2(testCase)
            actSolution = CX2(multicomplex([1,2,3,4,5,6,7,8]),2,3);
            expSolution = 7;
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testconsimulti(testCase)
            [actSolution1,actSolution2] = consimulti(multicomplex([1,2,3,4]),multicomplex([1,2,3,4,5,6,7,8]));
            expSolution = [multicomplex([1,2,3,4,0,0,0,0]),multicomplex([1,2,3,4,5,6,7,8])];
            testCase.verifyEqual([actSolution1,actSolution2],expSolution);
        end
        function testmodcheck(testCase)
            actSolution = modcheck(multicomplex([0.5,0.3,0.3,0.4]));
            expSolution = 'converge';
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testmodc(testCase)
            actSolution = modc(multicomplex([0.5,0.3,0.3,0.4]));
            expSolution = 0.739897714367898;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        function testyy(testCase)
            actSolution = yy(multicomplex([0.5,0.3,0.3,0.4]));
            expSolution = 0.768114574786861;
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',10^-10);
        end
        
        
    end
    
end 
