function [exact,a,b]=correct(i)
%In most of the 23 battery problems the integration bounds are
%0 and 1.
a=0;b=1;
switch i
  case 1
    exact= exp(1)-1;                  
  case 2
    exact= 0.7;                  
  case 3
    exact= 2/3;
  case 4                     
    exact=.47942822668880167;a=-1;
  case 5
    exact=1.5822329637296729331;a=-1;
  case 6
    exact=2/5;
  case 7           
    exact=2;
  case 8
    exact=.866972987339911038;   
  case 9
    exact= 1.15470053837925153;
  case 10    
    exact= log(2);
  case 11           
    exact= .379885493041722475;
  case 12
    exact= .777504634112248276;
  case 13
    exact= 0.00909863753916684291556;a=0.1;
  case 14
    exact= 0.5;b=10;
  case 15                    
    exact= 1;b=10;
  case 16                
    exact= .49936338107645674464;b=10;
  case 17
    exact= .11213930374163741027;a=0.01;
  case 18
    exact=.83867634269442961454;b=pi;
  case 19   
    exact=-.99999999999996446122;  
  case 20
    exact=1.5643964440690497731;a=-1;
  case 21 
    exact=.16349494301863722618; 
  case 22
    exact=-.63466518254339257343;  
  case 23
    exact=.013492485649467772692;
  otherwise
    exact=0;a=0;b=0; warning('Problem number is not between 1 and 23')
end
