classdef Polynomial 
   % Documentation example
   % A value class that implements a data type for polynonials
   
   properties(SetAccess = private,GetAccess = private)
      BernsCoeff
   end
   
   % Class methods
   methods(Access = public)
       
      % Constructor, cases: 1.- default constructor, 
      %                     2.- copy constructor, 
      %                     3.- control points constructor, 
      %                     4.- interpolator constructor,
      %                     5.- constructor from power representation
      function obj = Polynomial(varargin)
          if nargin==0
              obj.BernsCoeff = zeros(1); % Default constructor
          elseif nargin==1
              if isa(varargin{1},'Polynomial') % Copy constructor
                  polynomial = varargin{1};
                  obj.BernsCoeff = polynomial.BernsCoeff;
              else % Construct a Bernstein polynomial from its power representation
                  monomialCoeff = varargin{1};
                  obj.BernsCoeff = monomial2Bernstein(monomialCoeff);
              end
          elseif nargin==2
              if isa(varargin{1},'double') && isa(varargin{2},'double') % Construct a Bernstein polynomial by interpolation
                  parameters = varargin{1}; % interpolation nodes
                  interpolationPoints = varargin{2}; % interpolation points
                  obj.BernsCoeff = interpolation (parameters,interpolationPoints); 
              elseif isa(varargin{1},'double') && isa(varargin{2},'char')
                  coeff = varargin{1};
                  type = varargin{2};
                  switch type
                      case 'm' % Construct a Bernstein polynomial from its monomial representation
                          obj.BernsCoeff = monomial2Bernstein(coeff);
                      case 'c' % Construct a Bernstein polynomial from its control points
                          obj.BernsCoeff = coeff(:).';
                  end
              end
          end
      end % Polynomial
      
      % CHAR    CHAR(OBJ) creates a formated display of the polynom as
      % linear combination of Bernstein polynomials
      str = char(obj) 
      
      % DISP    DISP(OBJ) displays object obj in MATLAB syntax
      disp(obj)
             
      % GETDEGREE   GETDEGREE(OBJ) provides the apparent degree of the
      % polynomial (the true degree can be lower)
      n = getDegree(obj)
      
      % GETCOEFF    GETCOEFF(OBJ) provides the coefficients vector of the
      % polynomial 
      coeff = getCoeff(obj)

      % DOUBLE   DOUBLE(OBJ) returns the coefficients of the polynomial 
      % respect to the Bernstein basis
      coeff = double(obj)
      
      % DEGREEELEVATION     DEGREEELEVATION(OBJ,K) implements k degrees 
      % elevation of a polynomial represented in the Bernstein basis 
      poly = degreeElevation(obj,k)
      
      % DEGREEREDUCTION     DEGREEREDUCTION(OBJ) return a polynomial object
      % represented in the Bernstein basis of the true degree of the 
      % polynomial obj (it is no obvious the true degree of a polynomial 
      % represented in the Bernstein basis)
      poly = degreeReduction(obj)      

      % PLUS    PLUS(OBJ1,OBJ2) implements obj1+obj2 for Poylnomials
      r = plus(obj1,obj2)
      
      % MINUS   MINUS(obj1,obj2) implements obj1 - obj2 for Polynomials
      r = minus(obj1,obj2)
            
      % MTIMES   MTIMES(OBJ1,OBJ2) implements obj1 * obj2 for Polynomials
      r = mtimes(obj1,obj2)
            
      % MPOWER  MPOWER(OBJ,K) implements obj^k
      r = mpower(obj,k)

      % MRDIVIDE  MRDIVIDE(OBJ1,OBJ2) implements obj1/obj2
      [q,r] = mrdivide(obj1,obj2)                        
      
      % EVAL  EVAL(OBJ,X,PREC) evaluates obj at the points in x. varargin
      %       is an optional argument containing a pretended precision (if 
      %       this argument is not provided precision is taken as 1e-12)
      y = Eval(obj,x,varargin)

      % DIFF  DIFF(obj) provides the derivative of the polynomial obj like 
      %       an object of the class Polynomial
      r = diff(obj)
      
      % INTEGRAL  INTEGRAL(OBJ) provides the definite integral of the 
      %           polynomial obj in the interval [0,1]          
      q = integral(obj)
      
      % INTEGRATE  INTEGRATE(OBJ) provides the indefinite integral of the 
      %            polynomial obj like an object of the class Polynomial
      r = integrate(obj) 
      
      % PLOT  PLOT(OBJ,step) plots the Polynomial obj in the interval [0,1] 
      % from a uniform mesh points in this interval. If the step is not 
      % provided it is taken as 0.1
      plot(obj,varagin)         
    
      % SUBSREF Implementing the following syntax: 
      % obj([1 ...])
      % obj.BernsCoeff
      % obj.method
      % out = obj.method(args)
      % out = obj.method
      b = subsref(a,s)    
      
   end
   
end % classdef
