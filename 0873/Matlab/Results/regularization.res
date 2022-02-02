regularization

Problem: phillips.     Dimension: 300.    Delta: 2.999927e+000

Eigensolver: tcheigs_lstrs_gateway

--------------
init_up_bounds
--------------

alphaU = 2.7276e+002
deltaU = 2.3245e+001

--------
b_epairs
--------

alpha: 0.000000e+000

# converged eigenvalues: 1
lambda1: -6.913002e+001
nu1:     7.684500e-001
--------------
init_lo_bounds
--------------

alphaL = -9.6855e+001
deltaL = -6.9130e+001

----------
upd_deltaU
----------

not updating deltaU

   deltaN        deltaU
3.056199e+001 2.324526e+001

------------
adjust_alpha
------------

# of times alpha was adjusted: 0

LSTRS iteration: 0
------------
boundary_sol
------------

  | ||x|| - Delta|     eps_Delta * Delta
        2.1672              0.0300

||x|| = 0.8327  Delta = 2.9999

lambda1(alphak) = -69.1300

------------
interior_sol
------------

    ||u||       Delta*|nu|
6.399098e-001   2.3053

     lambda       -eps_Int
-6.913002e+001 -1.000000e-008

lambda > -eps_Int: 0
||x|| < Delta:     1

----------------
quasioptimal_sol
----------------
 
qopsol:       0

1. (l2-l1)*tau2^2*(Delta^2+1): 1.0000e+000, -2*eta*psi(xtilde): 0.0000e+000
2. (l2-l1)*tau2^2*(Delta^2+1): 1.0000e+000, -2*eta*psi(xtilde): 0.0000e+000
||x||: 8.327280e-001, lambda: -6.913002e+001
|||x||-Delta|/Delta:     7.224172e-001

------------
convergence
------------

iter:  0
boundary solution:        0
interior solution:        0
quasi-optimal solution:   0
small alpha interval:     0
max number of iterations: 0

-----------
inter_point
-----------

using 1st
 
--------------
upd_alpha_safe
--------------

alphaL = 0.0000e+000
alphaU = 2.7276e+002

-----------
safe_alpha1
-----------

not safeguarding alpha1
alpha = 2.519321e+002

--------
b_epairs
--------

alpha: 2.519321e+002

# converged eigenvalues: 2
lambda1: 7.791619e-005
nu1:     -2.460311e-008
lambda2: 3.066145e-004
nu2:     1.679031e-007

----------
upd_deltaU
----------

updating deltaU

   deltaN        deltaU
7.791619e-005 2.324526e+001

--------
b_epairs
--------

alpha: 1.259661e+002

# converged eigenvalues: 2
lambda1: -1.768532e+001
nu1:     4.969593e-001
lambda2: 6.161946e+000
nu2:     3.203367e-002

----------
upd_deltaU
----------

not updating deltaU

   deltaN        deltaU
2.942741e+001 7.791619e-005

------------
adjust_alpha
------------

# of times alpha was adjusted: 1

LSTRS iteration: 1
------------
boundary_sol
------------

  | ||x|| - Delta|     eps_Delta * Delta
        1.2538              0.0300

||x|| = 1.7462  Delta = 2.9999

lambda1(alphak) = -17.6853

------------
interior_sol
------------

    ||u||       Delta*|nu|
8.677739e-001   1.4908

     lambda       -eps_Int
-1.768532e+001 -1.000000e-008

lambda > -eps_Int: 0
||x|| < Delta:     1

----------------
quasioptimal_sol
----------------
 
tau1 = 5.8401e-001
tau2 = 8.1174e-001

tau1 = 6.8340e-001
tau2 = -7.3005e-001

nu1^2 + nu2^2: 2.4799e-001 1/(1+Delta^2): 1.0000e-001

product >= 1: 1

qopsol:       0

1. (l2-l1)*tau2^2*(Delta^2+1): 1.5713e+002, -2*eta*psi(xtilde): 1.4568e-006
2. (l2-l1)*tau2^2*(Delta^2+1): 1.2709e+002, -2*eta*psi(xtilde): 1.7572e-006
||x||: 1.746167e+000, lambda: -1.768532e+001
|||x||-Delta|/Delta:     4.179302e-001

------------
convergence
------------

iter:  1
boundary solution:        0
interior solution:        0
quasi-optimal solution:   0
small alpha interval:     0
max number of iterations: 0

-----------
inter_point
-----------

using 1st
 
--------------
upd_alpha_safe
--------------

alphaL = 1.2597e+002
alphaU = 2.5193e+002

---------
interpol2
---------
 
safeguarding lambdabar, greater than deltaU

-----------
safe_alphak
-----------

not safeguarding alphak
alpha = 2.301234e+002

--------
b_epairs
--------

alpha: 2.301234e+002

# converged eigenvalues: 2
lambda1: -3.757879e-001
nu1:     3.224773e-001
lambda2: 2.730958e-005
nu2:     -1.550805e-007

----------
upd_deltaU
----------

not updating deltaU

   deltaN        deltaU
2.637618e+001 7.791619e-005

------------
adjust_alpha
------------

# of times alpha was adjusted: 0

LSTRS iteration: 2
------------
boundary_sol
------------

  | ||x|| - Delta|     eps_Delta * Delta
        0.0646              0.0300

||x|| = 2.9353  Delta = 2.9999

lambda1(alphak) = -0.3758

------------
interior_sol
------------

    ||u||       Delta*|nu|
9.465772e-001   0.9674

     lambda       -eps_Int
-3.757879e-001 -1.000000e-008

lambda > -eps_Int: 0
||x|| < Delta:     1

----------------
quasioptimal_sol
----------------
 
tau1 = 9.8064e-001
tau2 = 1.9581e-001

tau1 = 9.8064e-001
tau2 = -1.9581e-001

nu1^2 + nu2^2: 1.0399e-001 1/(1+Delta^2): 1.0000e-001

product >= 1: 1

qopsol:       0

1. (l2-l1)*tau2^2*(Delta^2+1): 1.4409e-001, -2*eta*psi(xtilde): 2.3374e-006
2. (l2-l1)*tau2^2*(Delta^2+1): 1.4409e-001, -2*eta*psi(xtilde): 2.3374e-006
||x||: 2.935330e+000, lambda: -3.757879e-001
|||x||-Delta|/Delta:     2.153292e-002

------------
convergence
------------

iter:  2
boundary solution:        0
interior solution:        0
quasi-optimal solution:   0
small alpha interval:     0
max number of iterations: 0

-----------
inter_point
-----------

using 1st
 
--------------
upd_alpha_safe
--------------

alphaL = 2.3012e+002
alphaU = 2.5193e+002

---------
interpol2
---------
 
safeguarding lambdabar, greater than deltaU

-----------
safe_alphak
-----------

not safeguarding alphak
alpha = 2.337458e+002

--------
b_epairs
--------

alpha: 2.337458e+002

# converged eigenvalues: 2
lambda1: -5.952830e-003
nu1:     -3.163805e-001
lambda2: 7.790953e-005
nu2:     -7.398133e-006

----------
upd_deltaU
----------

not updating deltaU

   deltaN        deltaU
2.599435e+001 7.791619e-005

------------
adjust_alpha
------------

# of times alpha was adjusted: 0

LSTRS iteration: 3
------------
boundary_sol
------------

  | ||x|| - Delta|     eps_Delta * Delta
        0.0015              0.0300

||x|| = 2.9984  Delta = 2.9999

lambda1(alphak) = -0.0060

------------
interior_sol
------------

    ||u||       Delta*|nu|
9.486324e-001   0.9491

     lambda       -eps_Int
-5.952830e-003 -1.000000e-008

lambda > -eps_Int: 0
||x|| < Delta:     1

----------------
quasioptimal_sol
----------------
 
tau1 = -9.9954e-001
tau2 = -3.0381e-002

tau1 = -9.9954e-001
tau2 = 3.0334e-002

nu1^2 + nu2^2: 1.0010e-001 1/(1+Delta^2): 1.0000e-001

product >= 1: 1

qopsol:       0

1. (l2-l1)*tau2^2*(Delta^2+1): 5.5660e-005, -2*eta*psi(xtilde): 2.3381e-006
2. (l2-l1)*tau2^2*(Delta^2+1): 5.5489e-005, -2*eta*psi(xtilde): 2.3381e-006
||x||: 2.998391e+000, lambda: -5.952830e-003
|||x||-Delta|/Delta:     5.121129e-004

------------
convergence
------------

iter:  3
boundary solution:        1
interior solution:        0
quasi-optimal solution:   0
small alpha interval:     0
max number of iterations: 0

-------
output
-------
lambda1: -5.952830e-003
nu1:     -3.163805e-001
lambda2: 7.790953e-005
nu2:     -7.398133e-006


Number of LSTRS Iterations:     4
Number of calls to eigensolver: 5
Number of MV products:          721


(||x||-Delta)/Delta: 5.121129e-004

lambda: -5.952830e-003

||g + (H-lambda* I)x||/||g|| = 5.639093e-015

The vector x is a Boundary Solution
