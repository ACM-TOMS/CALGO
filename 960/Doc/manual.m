% This test file contains all the examples in the Manual.pdf file in the Doc 
% directory

clc
clear

disp(' ');
disp('CONSTRUCTING OBJECTS OF THE CLASS');
disp(' ');
disp('1.-The default constructor creates the zero degree polynomial whose');
disp('value is zero at any point');
disp(' ');
disp('>> poly1 = Polynomial()');

poly1 = Polynomial()

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');


disp('2.- If p(x)=a0+a1x+...+anx^n, by using Polynomial([a0,a1,...,an])');
disp('or Polynomial([a0,a1,...,an],''m''), we contruct this polynomial');
disp('represented in the basis formed by the Bernstein polynomials of'); 
disp('degree n. So, let us construct the polynomial p(t)=1+2x+3x^2 in both');
disp('ways');
disp(' ');
disp('>> poly2=Polynomial([1,2,3])');

poly2 = Polynomial([1,2,3])

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp('>> poly2=Polynomial([1,2,3],''m'')');

poly2 = Polynomial([1,2,3],'m')

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('3.- Given the interpolation conditions p(xi) = qi, i = 0,...,n,');
disp('with 0<x0<x1<...<xn<1, we can construct the interpolating polynomial'); 
disp('represented in the basis formed by the Bernstein polynomials of degree'); 
disp('n by invoking Polynomial([x0,x1,...,xn],[q0,q1,...,qn]). Let us'); 
disp('cosntruct the degree 2 polynomial satisfying p(0.25)=1, p(0.5)=-2'); 
disp('and p(0.75)=-3:');
disp(' ');
disp('>> poly3=Polynomial([0.25,0.5,0.75],[1,-2,3])');

poly3 = Polynomial([0.25,0.5,0.75],[1,-2,3])

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('4.- If the coefficients ci of the polynomial with respect to the basis');
disp('formed by the Bernstein polynomials of degree n are known, then the');
disp('polynomial can be constructed by using Polynomial([c0,c1,...,cn],''c'').'); 
disp(' So, let us construct the polynomial with Bernstein coefficients 0,');
disp('1,2 and 3:');
disp(' ');
disp('>> poly4 = Polynomial(0:1:3,''c'')');

poly4 = Polynomial(0:1:3,'c')

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('5.- Finally, the copy constructor, which can be invoked by performing');
disp('an assignement');
disp(' ');
disp('>> poly5 = poly2');

poly5 = poly2

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');

disp('or, equivalently, by using explicitly the constructor of the class');
disp(' ');
disp('>> poly5 = Polynomial(poly2)');

poly5 = Polynomial(poly2)

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');
disp(' ');
disp(' ');



disp('USING THE METHODS OF THE CLASS');
disp(' ');
disp('1. The method getDegree returns the apparent degree of the polynomial');
disp(' ');
disp('>> poly5.getDegree()');

poly5.getDegree()

disp('>> poly4.getDegree()');

poly4.getDegree()

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');


disp('2. The methods getCoeff and double return the coefficients of the'); 
disp('polynomial respect to the Bernstein basis');
disp(' ');
disp('>> poly5.getCoeff()');

poly5.getCoeff()

disp('>> poly4.double()');

poly4.double()

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('3. A polynomial represented in a Bernstein basis of a certain degree');
disp('can be also represented exactly in a Bernstein basis of a greater');
disp('degree. The method degreeElevation performs this degree elevation.');
disp('As example let us express the degree 2 polynomial stored in poly5');
disp('using the Bernstein basis of degree 4(=2+2):');
disp(' ');
disp('>> poly5 = poly5.degreeElevation(2)');

poly5 = poly5.degreeElevation(2)

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');

disp('4. Although a polynomial is represented in a Bernstein basis of a ');
disp('certain degree n, the true degree of the polynomial can be lower. ');
disp('The method degreeReduction checks the true degree of a polynomial'); 
disp('and returns the same polynomial represented in the Bernstein basis');
disp('of the lowest possible degree. Taking into account that now poly5');
disp('is a polynomial of degree at most 2 represented in a Bernstein basis');
disp('of degree 4 let us reduce its degree:');
disp(' ');
disp('>> poly5 = poly5.degreeReduction()');

poly5 = poly5.degreeReduction()

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('5. The usual arithmetic operations can be performed with the usual');
disp('operators + (sum), - (substraction), * (product), / (division) and');
disp('^ (power). Let us see some examples:');
disp(' ');
disp('>> poly2');

poly2

disp('>> poly4');

poly4

disp('>> poly2+poly4');

poly2+poly4

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');

disp('>> poly2-poly4');

poly2-poly4

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');

disp('>> poly2*poly4');

poly2*poly4

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');

disp('>> [q,r] = poly2/poly4');

[q,r] = poly2/poly4

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');

disp('>> poly6 = poly4*q+r');

poly6 = poly4*q+r

disp('>> poly6.degreeReduction()');

poly6.degreeReduction()

disp('>> poly2');

poly2

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');

disp('>> poly2^2');

poly2^2

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('6. The method Eval evaluates a polynomial. It can be invoked in two');
disp('different ways: Eval(x) and Eval(x,prec). The first one evaluates'); 
disp('the polynomial at the points in x with a default pretended percision'); 
disp('of 1e-12, whereas the second one performs the same evaluation with');
disp('a pretended precision given by prec. Let us see some examples:');
disp(' ');
disp('>> format long e');

format short e

disp('>> poly4.Eval([0:0.25:1])');

poly4.Eval([0:0.25:1])

disp(' ');
disp('As we can observe, the method returns for the row vector of 5 points'); 
disp('a 3x5 matrix, where the first row consists of the evaluations of the');
disp('polynomial at the corresponding points in x with a pretended precision');
disp('of 1e-12, the second row provides, when possible, upper bounds of the');
disp('relative error for the evaluations, and the third row consists of');
disp('flags (a 1 flag means that relative error is lower than prec, a 0 ');
disp('flag means that either relative error is not less than prec or it is'); 
disp('not known about.'); 
disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp('If the vector of points x is a column vector instead of a row vector');
disp('we obtain the traspose matrix:');

disp('>> poly4.Eval([0;0.25;0.5;0.75;1])');

poly4.Eval([0;0.25;0.5;0.75;1])

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp('In case that we are not satisfied with the precision obtained we can');
disp('require a greater (or lower) precision by using the second argument');
disp('prec of the method:')
disp(' ');

disp('>> poly4.Eval([0:0.25:1],5*1e-16)');

poly4.Eval([0:0.25:1],5*1e-16)

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');

disp('Evaluations of the polynomials can also be obtained with obj(x) or obj(x,prec):');
disp(' ');
disp('>> poly4(0:0.25:1)');

poly4(0:0.25:1)

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('7. The method diff returns the derivate of a polynomial. For example,');
disp(' ');

disp('>> poly4');

poly4

disp(' ');
disp('>> poly7 = poly4.diff()');

poly7 = poly4.diff()

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('8. The method integrate returns the indefinite integral of a polynomial.')
disp(' ');

disp('>> poly7.integrate()');

poly7.integrate()

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('9. integral returns the definite integral of a polynomial in the'); 
disp('interval [0,1]');
disp(' ');
disp('>> poly7.integral()');

poly7.integral()

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');
disp(' ');



disp('10. The method plot display a graphic of the polynomial in [0,1].');
disp('This method can be invoked in two ways: plot() or plot(step). The');
disp('first one, makes the graphic taking a uniform partition of the ');
disp('interval [0,1] with step 0.01, whereas the second one takes a uniform');
disp('partition according to step. Let us see a couple of examples:');
disp(' ');

disp('>> figure');

figure;

disp('>> poly5.plot()');

poly5.plot()

disp(' ');
disp('Press any key to continue');
pause;
disp(' ');

disp('>> figure');

figure;

disp('>> poly5.plot(0.25)');

poly5.plot(0.25)
