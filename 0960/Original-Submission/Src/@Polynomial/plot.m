% PLOT  PLOT(OBJ,step) plots the Polynomial obj in the interval [0,1] from
% a uniform mesh points in this interval. If the step is not provided it is
% taken as 0.1
function plot(obj,varargin)
   if nargin==1
       step = 0.01;
   elseif nargin==2
       step = varargin{1};
       if step>0.5
           error('Polynomial.plot: step must be lower or equal to 0.5');
       end
   else
       error('Polynomialplot: too input argumentes');
   end
   x = 0:step:1; 
   y = Eval(obj,x);
   plot(x,y(1,:));
   c = char(obj);
   title(['y = ' c])
   xlabel('X')
   ylabel('Y','Rotation',0)
   dist = (max(y(1,:))-min(y(1,:)))*0.05;
   axis([-0.05,1.05,min(y(1,:))-dist,max(y(1,:))+dist]);
   grid on
end 
