with Floating_Point_Numbers;   use Floating_Point_Numbers;
with Mathematical_Functions;   use Mathematical_Functions;

package body Complex_Numbers is

  TWO : constant double_float := double_float(2);

-- real => complex :

  function CMPLX ( X : double_float ) return double_complex is
  begin
    return ( X, ZERO );
  end CMPLX;

-- real * real => complex :

  function CMPLX ( X,Y : double_float ) return double_complex is
  begin
    return ( X, Y );
  end CMPLX;

-- real * argument => complex :

  function CMPLX_POLAR ( R : double_float; A : double_float )
                       return double_complex is 
  begin
    return ( R*SIN(A), R*COS(A) );
  end CMPLX_POLAR;

-- complex => real :

  function RADIUS ( X, Y : double_float) return double_float is
  begin
    return Sqrt(x**2+y**2);
  end RADIUS;

  function "ABS" ( X : double_complex ) return double_float is 
  begin 
    return RADIUS(X.RE,X.IM);
  end "ABS";

  function REAL_PART ( X : double_complex ) return double_float is
  begin
    return X.RE;
  end REAL_PART;

  function IMAG_PART ( X : double_complex ) return double_float is
  begin
    return X.IM;
  end IMAG_PART;

-- complex => argument :

  function ANGLE ( X, Y : double_float ) return double_float is
  begin
    return Mathematical_Functions.ANGLE(x,y);
  end ANGLE;
 
  function ARGUMENT ( X : double_complex ) return double_float is
  begin
    return ANGLE(X.IM,X.RE);
  end ARGUMENT;

-- complex * real => complex :

  function "+"  ( X : double_complex; Y : double_float )
                return double_complex is
  begin
    return (X.RE + Y, X.IM);
  end "+";

  function "-" ( X : double_complex; Y : double_float )
                return double_complex is
  begin
    return (X.RE - Y, X.IM);
  end "-";

  function "*"  ( X : double_complex; Y : double_float )
                return double_complex is
  begin
    return ( X.RE*Y, X.IM*Y );
  end "*";

  function "/"  ( X : double_complex; Y : double_float )
                return double_complex is
  begin
    return ( X.RE/Y, X.IM/Y );
  end "/";

-- real * complex => complex :

  function "+"  ( X : double_float; Y : double_complex )
                return double_complex is
  begin
    return ( X + Y.RE, +Y.IM );
  end "+";

  function "-"  ( X : double_float; Y : double_complex )
                return double_complex is
  begin
    return ( X - Y.RE, -Y.IM );
  end "-";

  function "*"  ( X : double_float; Y : double_complex )
                return double_complex is
  begin
    return ( X*Y.RE, X*Y.IM );
  end "*";

  function "/"  ( X : double_float; Y : double_complex )
                return double_complex is

    R: double_float renames X;
    U: double_float renames Y.RE;
    T: double_float renames Y.IM;

  begin
    if      U = ZERO then
      return ( R / T ,
               ZERO   ) ;
    elsif   T = ZERO then
      return (  ZERO  ,
               -R / U  ) ;
    elsif  abs(T) < abs(U) then
      return ( ( R*(T/U) ) / ( T*(T/U) + U ) ,
               ( - R     ) / ( T*(T/U) + U )  ) ;
    elsif  abs(T) > abs(U) then
      return ( (   R       ) / ( T + U*(U/T) ) ,
               ( - R*(U/T) ) / ( T + U*(U/T) )  ) ;
    elsif   T = U  then
      return ( (   R ) / ( TWO*T ) , 
               ( - R ) / ( TWO*T )  ) ;
    else -- T = -U then
      return ( ( - R ) / ( TWO*T ) , 
               ( - R ) / ( TWO*T )  ) ;
    end if;
  end "/";

-- complex => complex :

  function "+"   ( X : double_complex ) return double_complex is
  begin
    return X; -- null function
  end "+";

  function "-"   ( X : double_complex ) return double_complex is
  begin
    return ( -X.RE, -X.IM );
  end "-";

-- complex * complex => complex :

  function "+"  ( X,Y : double_complex ) return double_complex is
  begin
    return ( X.RE + Y.RE, X.IM + Y.IM );
  end "+";

  function "-"  ( X,Y : double_complex ) return double_complex is
  begin
    return ( X.RE - Y.RE, X.IM -Y.IM );
  end "-";

  function "*"  ( X,Y : double_complex ) return double_complex is
    R: double_float renames X.RE;
    S: double_float renames X.IM;
    T: double_float renames Y.RE;
    U: double_float renames Y.IM;
  begin
    if      U = ZERO then
      return ( R * T ,
               S * T  ) ;
    elsif   T = ZERO then
      return ( -S * U ,
                R * U  ) ;
    elsif  abs(T) < abs(U) then
      return ( ( R*(T/U) - S ) * U ,
               ( S*(T/U) + R ) * U ) ;
    elsif  abs(T) > abs(U) then
      return ( ( R - S*(U/T) ) * T ,
               ( S + R*(U/T) ) * T ) ;
    elsif   T = U  then
      return ( ( R - S ) * T , 
               ( S + R ) * T  ) ;
    else -- T = -U then
      return ( ( R + S ) * T , 
               ( S - R ) * T  ) ;
    end if;
  end "*";

  function "/"  ( X,Y : double_complex ) return double_complex is

    R: double_float renames X.RE;
    S: double_float renames X.IM;
    T: double_float renames Y.RE;
    U: double_float renames Y.IM;

  begin
    if      U = ZERO then
      return ( R / T ,
               S / T  ) ;
    elsif   T = ZERO then
      return (  S / U , -R / U  ) ;
    elsif  abs(T) < abs(U) then
      return ( ( R*(T/U) + S ) / ( T*(T/U) + U ) ,
               ( S*(T/U) - R ) / ( T*(T/U) + U )  ) ;
    elsif  abs(T) > abs(U) then
      return ( ( R + S*(U/T) ) / ( T + U*(U/T) ) ,
               ( S - R*(U/T) ) / ( T + U*(U/T) )  ) ;
    elsif   T = U  then
      return ( ( R + S ) / ( TWO*T ) , 
               ( S - R ) / ( TWO*T )  ) ;
    else -- T = -U then
      return ( ( R - S ) / ( TWO*T ) , 
               ( S + R ) / ( TWO*T )  ) ;
    end if;
  end "/";

  function "**" ( c : double_complex; m : integer ) return double_complex is
    r : double_complex;
  begin
    if m = 0
     then return CMPLX(1.0);
     elsif m > 0
         then r := c;
              for j in 2..m loop
                r := c*r;
              end loop;
         else r := CMPLX(1.0);
              for j in 1..(-m) loop
                r := r/c;
              end loop;
    end if;
    return r;
  end "**";

  function de_Moivre (n,j : natural; c : double_complex)
                     return double_complex is

    arg,radius_c,angle_c : double_float;
    tmp : double_complex;

  begin
    arg := (2.0 * PI * double_float(j)) / double_float(n);
    radius_c := RADIUS(c)**(1.0/double_float(n));
    angle_c := ANGLE(c)/double_float(n);
    tmp := CMPLX(radius_c)*CMPLX(COS(angle_c),SIN(angle_c));
    return CMPLX(COS(arg),SIN(arg))*tmp;
  end de_Moivre;

end Complex_Numbers;
