* example1.gms
* This file contains the GAMS scalar format description of the problem 
*
* minimize objvar = -2*x1 +3*x2 -2*x3
* subject to 
*       x1^2 + 3*x2^2 -2*x2*x3 +3*x3^2 -17*x1 +8*x2 -14*x3 >= -19, 
*       x1 + 2*x2 + x3 <= 5, 
*       5*x2 + 2*x3 = 7
*       0 <= x1 <= 2, 0 <= x2 <= 1. 
*
* To solve this problem by sparsePOP.m:
* >> param.relaxOrder = 3;
* >> sparsePOP('example1.gms',param); 
* 
* This problem is also described in terms of the SparsePOP format 
* in the file example1.m.  See Section 3 of the manual.
*
* To obtain a tight bound for the optimal objective value by the function
* sparsePOP.m, set the parameter param.relaxOrder = 3.

* The description consists of 5 parts except comment lines
* starting the character '*'. The 5 parts are:
* < List of the names of variables >
* < List of the names of nonnegative variables >
* < List of the names of constraints >
* < The description of constraints >
* < Lower and upper bounds of variables  >

* < List of the names of variables >
Variables  x1,x2,x3,objvar;
* 'objvar' represents the value of the objective function.

* < List of the names of nonnegative variables >
Positive Variables x1, x2;

* < List of the names of constraints >
Equations  e1,e2,e3,e4;

* < The description of constraints >
* Each line should start with the name of a constraint in the list of names
* of constraints,  followed by '.. '. The symbols '*', '+', '-', '^', '=G='
* (not less than), '=E=' (equal to) and '=L=' (not larger than) can be used 
* in addition to the variables in the list of the names of variables and real 
* numbers. The right-hand side of an inequality or an equality should be a 
* single constant. One constraint can be described in more than one lines; 
* for example, 
* e2..    - 17*x1 + 8*x2 - 14*x3 +6*x1^2 + 3*x2^2 - 2*x2*x3 + 3*x3^2 =G= -19;
* is equivalent to
* e2..     - 17*x1 + 8*x2 - 14*x3 +6*x1^2
*                + 3*x2^2 - 2*x2*x3 + 3*x3^2 =G= -19;
* Note that the first letter of a line can not be '*' except comment lines.

* minimize objvar = -2*x1 +3*x2 -2*x3
e1..    2*x1 - 3*x2 + 2*x3 + objvar =E= 0;

* 6*x1^2 + 3*x2^2 -2*x2*x3 +3*x3^2 -17*x1 +8*x2 -14*x3 >= -19
e2..    - 17*x1 + 8*x2 - 14*x3 +6*x1^2 + 3*x2^2 - 2*x2*x3 + 3*x3^2 =G= -19;

* x1 + 2*x2 + x3 <= 5
e3..    x1 + 2*x2 + x3 =L= 5;

* 5*x2 + 2*x3 = 7
e4..    5*x2 + 2*x3 =E= 7;

* < Lower and upper bounds on variables  >
* Each line should contain exactly one bound; 
* For 0.5 <= x3 <= 2, we set 
* x3.lo = 0.5; 
* x3.up = 2; 
* A line such that 'x3.lo = 0.5; x3.up = 2;' is not allowed.

* x1 <= 2
x1.up = 2;

* x2 <= 1
x2.up = 1;

* end of example1.gms
