 function options = disodeset( name1, value1,name2,value2,name3,value3,name4,value4,name5,value5,name6,value6,name7,value7)
%  disodeset Create/alter DISOD45E OPTIONS structure.
%     OPTIONS = disodeset('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%     options structure OPTIONS in which the named properties have the
%     specified values. Any unspecified properties have default values.
%     
%     
%  disodeset PROPERTIES
%     
%  RelTol - Relative error tolerance  [ positive scalar {1e-3} ]
%     This scalar applies to all components of the solution vector, and
%     defaults to 1e-3 (0.1% accuracy).  The estimated error in
%     each integration step satisfies ||est|| <= RelTol*||y||+AbsTol.
%  
%  AbsTol - Absolute error tolerance  [ positive scalar or vector {1e-6} ]
%     AbsTol defaults to 1e-6 in all solvers. See RelTol.
%     
%  Refine - Output refinement factor  [ positive integer ]
%     This property increases the number of output points by the specified
%     factor producing smoother output. Refine defaults to 4.
%     
%  ActionSwitch - Installable output function  [ function_handle ]
%     This output function is called by the solver once a discontinuity
%     has been found, if the corresponding value of the parameter "isterminal"
%     is negative. ActionSwitch defaults to [].
%  
%  InitialStep - Suggested initial step size  [ positive scalar ]
%     The solver will try this first.  By default the solvers determine an
%     initial step size automatically. 
%
%  Gradient -  Function handle for the gradien of the event function
%
%  EventControl - Type of control for the detection of discontinuities
%     0   Existence of discontinuity is checked at every step and every
%         stage of failed steps
%     1   Existence of discontinuity is checked at every stage of every
%           1  step.
%     k   Existence of discontinuity is checked at every stage of every
%         step and at 6*k uniformly distributed points inside every step.
%
options=struct('AbsTol',[],'RelTol',[],'Gradient',[],'InitialStep',[],'EventControl',[],'ActionSwitch',[],'Refine',[]);
if nargin==2
    options=setfield(options,name1,value1);
elseif nargin==4
    options=setfield(options,name1,value1);
    options=setfield(options,name2,value2);
elseif nargin ==6
    options=setfield(options,name1,value1);
    options=setfield(options,name2,value2);
    options=setfield(options,name3,value3);
elseif nargin ==8
    options=setfield(options,name1,value1);
    options=setfield(options,name2,value2);
    options=setfield(options,name3,value3);
    options=setfield(options,name4,value4);
elseif nargin ==10
    options=setfield(options,name1,value1);
    options=setfield(options,name2,value2);
    options=setfield(options,name3,value3);
    options=setfield(options,name4,value4);
    options=setfield(options,name5,value5);
elseif nargin ==12
    options=setfield(options,name1,value1);
    options=setfield(options,name2,value2);
    options=setfield(options,name3,value3);
    options=setfield(options,name4,value4);
    options=setfield(options,name5,value5);
    options=setfield(options,name6,value6);
elseif nargin >=14
    options=setfield(options,name1,value1);
    options=setfield(options,name2,value2);
    options=setfield(options,name3,value3);
    options=setfield(options,name4,value4);
    options=setfield(options,name5,value5);
    options=setfield(options,name6,value6);
    options=setfield(options,name7,value7);

end

