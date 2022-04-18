#When a calculation fails, we will increase Digits by d_step
d_step := 10 :

S := proc( n :: integer , x :: numeric , d :: integer )
#Calculates S_n(x) to at least d decimal digits using formula (28) in
#the paper. We only need to make these calculations once, so performance 
#doesn't matter. We will use brute force to get the correct values, 
#simply increasing the number of digits until the answer converges.

  local j   , p   :
  local pos , neg :
  local t   , tmp :
  local f         :
  local S         :

  local results := Array( 0..1 ) :

  #-----

  #This will probably be sufficient ...
  Digits := 2 * d :

  if ( x = 0 ) then
    return evalf( 1 / ( n + 1/2 ) , d ) , Digits :
  end if :
 
  #... but this loop repeats until the required accuracy is achieved, just in case.
  for p from 0 by 1 do 

    t := 1 / evalf( n + 1/2 ) :
    
    if ( n > 0 ) then 
      pos := t :
      neg := 0 :
    else
      pos := t :
      neg := 0 :
    end if :

    for j from 1 by 1 do 

      #Calculate term
      f := x / evalf( n + j + 1/2 ) :
      t := -t * f :

      if ( t > 0 ) then
 
        tmp := pos + t :

        #If including the current term has no effect we will stop the
        #loop.  Note that smaller terms may still affect the negative
        #part of the sum, but these would be lost when we combine the
        #two parts at the end.

        if ( tmp = pos ) then
          if ( j >= -n and abs( f ) <= 1 ) then
            break : 
          end if :
        else
          pos := tmp :
        end if

      else

        tmp := neg + t :

        #If including the current term has no effect we will stop the
        #loop.  Note that smaller terms may still affect the positive
        #part of the sum, but these would be lost when we combine the
        #two parts at the end.

        if ( tmp = neg ) then
          if ( j >= -n and abs( f ) <= 1 ) then
            break : 
          end if :
        else
          neg := tmp :
        end if

      end if :

    end do :

    #-----

    #Note: this addition must be not be performed inside an 
    #evalf() statement with a reduced number of digits or else
    #rounding will take place too early.
    S := pos + neg : 

    #Unlikely --- all digits lost due to cancellation. Try again at higher 
    #precision. This means that the function cannot return 0 as a result, 
    #so we are assuming that the locations of the zeros are irrational.
    if ( S = 0 ) then
      Digits := 2 * Digits :
      repeat :
    end if :

    if ( p <= 1 ) then

      results[p] := evalf( S , d ) : #OK to round the result now.

    else

      results[p mod 2] := evalf( S , d ) :

      if ( results[0] = results[1] ) then
        break :
      end if :

    end if :

    #-----

    Digits := Digits + d_step :

    #-----

  end do :

  #-----

  return evalf( results[0] , d ) , Digits :

#-----

end proc :    

#-----

xS := proc( n :: integer , x :: numeric , d :: integer ) 
#Calculates exp( x ) * S_n(x) to at least d decimal digits using formula 
#(12) in the paper. We only need to make these calculations once, so 
#performance doesn't matter. We will use brute force to get the correct 
#values, simply increasing the number of digits until the answer converges.

  local j   , p     : 
  local xjj , f , t :
  local pls , mns   :
  local tmp , S     :

  local results := Array( 0..1 ) :

  #-----

  #This will probably be sufficient ...
  Digits := 2 * d :

  if ( x = 0 ) then
    return evalf( 1 / ( n + 1/2 ) , d ) , Digits :
  end if :
 
  #... but this loop repeats until the required accuracy is achieved, just in case.
  for p from 1 by 1 do    

    xjj := 1 :  # x^j / j!

    if ( n >= 0 ) then
      pls := evalf( 1 / ( n + 1/2 ) )  :
      mns := 0 :
    else
      pls := 0 :
      mns := evalf( 1 / ( n + 1/2 ) ) :
    end if :

    #-----

    for j from 1 by 1 do

      xjj := xjj * x / j :
        t := evalf( xjj / ( j + n + 1/2 ) ) :

      if ( t > 0 ) then

        tmp := pls + t :

        if ( tmp = pls ) then

          #Check that terms subsequent terms are smaller
          if ( j >= max( -n , abs( x ) ) ) then   
            break :
          end if :

        else
          pls := tmp :
        end if

      else

        tmp := mns + t :
  
        if ( tmp = mns ) then

          #Check that terms subsequent terms are smaller
          if ( j >= max( -n , abs( x ) ) ) then   
            break :
          end if :

        else
          mns := mns + t :
        end if :

      end if :

    end do :

    #-----

    #Note: this addition must be not be performed inside an 
    #evalf() statement with a reduced number of digits or else
    #rounding will take place too early.
    S := pls + mns : 

    #Unlikely --- all digits lost due to cancellation. Try again at higher 
    #precision. This means that the function cannot return 0 as a result, 
    #so we are assuming that the locations of the zeros are irrational.
    if ( S = 0 ) then
      Digits := 2 * Digits :
      repeat :
    end if :      

    if ( p <= 1 ) then

      results[p] := evalf( S , d ) : #OK to round the result now.

    else

      results[p mod 2] := evalf( S , d ) :

      if ( results[0] = results[1] ) then
        break :
      end if :

    end if :

    #-----

    Digits := Digits + d_step :

    #-----

  end do :

  #-----

  return evalf( results[0] , d ) , Digits :

  #-----

end proc :

#-----

double_calc_S := proc( n :: integer , x :: numeric , d :: integer ) 
#Computes Sn( x ) to using two different methods, increasing
#precision until the results agree to d significant digits.

  local dl := max( d , 1 ) :

  local s1 , s2 :

  #-----

  do

    Digits := dl :

    s1 := S( n , x , dl )[1] :

    #One extra digit needed to negate rounding error in multiplication.
    s2 := evalf( exp( -x ) , dl + 1 ) * xS( n , x , dl + 1 )[1] :

    if ( evalf( s1 , d ) = evalf( s2 , d ) ) then
      return evalf( s1 , d ) :
    else 
      dl := dl + d_step :
    end if :

  end do :

end proc :

#-----

S_root := proc( n :: integer , d :: integer )
#Calculates x such that | x - xn | < xn * 10**(-d), where x = xn 
#is the zero of S( n , x ), and n < 0. Uses S if scale is false,
#and xS if scale is true.

  local x , xmin , xmax :
  local sval :
  local eps  :

  #-----

  if ( n >= 0 ) then
    error "S_root: invalid n value (must be negative)" :
  end if :

  #-----

  #If the value of d is very small, setting Digits := d could
  #lead to the initial values of xmin and xmax being the same.
  Digits :=  max( d , ceil( log10( abs( n ) ) ) ) + 2 :

  #Initial guess based on asymptotics
  x      := -evalf( n + 1/6 ) : 
  xmin   := x - 0.1 :
  xmax   := x + 0.1 :
  
  #-----

  #Make sure that S( n , xmin ) < 0 and S( n , xmax ) > 0.

  #We only need to know whether S is positive or negative;
  #therefore we don't require much precision here.

  do

    sval := double_calc_S( n , xmin , 4 ) :

    if ( sval < 0 ) then
      break :
    end if :

    xmin := xmin - 0.1 :

    if ( xmin < 0 ) then
      xmin = 0 :
      break    :
    end if :

  end do :

  #-----

  do

    sval := double_calc_S( n , xmax , 4 ) :

    if ( sval > 0 ) then
      break :
    end if :

    xmax := xmax + 0.1 :

  end do :

  #-----

  #Refine xmin and xmax using a bisection search. We want 
  #| x - xn | / xn <= 0.5 * 10**(-d) where xn is the exact
  #location of the root. A sufficient condition for this is
  #( xmax - xmin ) / xmin <= 10**(-d). Now the value of 
  #( xmax - xmin ) after j steps is known, but xmin itself 
  #is not. Therefore we won't bother trying to predict the 
  #number of steps needed.
  eps := 10**(-d) :

  do 

    x    := ( xmin + xmax ) / 2   : 
    sval := double_calc_S( n , x , 4 ) :

    #Note: the function S cannot return 0.
    if ( sval > 0 ) then
      xmax := x :
    else
      xmin := x :
    end if :

    if ( xmax / xmin - 1 < eps ) then
      #Search is complete; minimize maximum
      #error by moving to interval centre.
      x := ( xmax + xmin ) / 2 : 
      break :
    end if :

  end do :

  #-----

  return evalf( x , d ) :

  #-----

end proc :

#-----

LibraryTools[Save]( 'd_step', 'S' , 'xS' , 'double_calc_S' , 'S_root' , "./" ) :

#-----

#End of file 'in_gamma.mpl'`