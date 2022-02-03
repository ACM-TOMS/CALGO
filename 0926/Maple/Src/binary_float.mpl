#Provides mechanisms for converting Maple's decimal
#floats to and from binary floats of arbitrary length,
#and some simple utilities for use with binary floats.

#Binary floats are represented as an ordered triple of 
#integers ( sg , mantis , xp ), where sg is +1 or -1.

#-----

#These parameters define the binary float
nbits         := 64 :
mantis_length := 53 : # Includes hidden bit; this is p in IEEE754
emax          := 1023 :

#These are calculated from the above.
emin          := 1 - emax :
bias          := emax :
mantis_min    := 2^(mantis_length-1) :
mantis_max    := 2 * mantis_min - 1 :

#Number of decimal significant figures that
#is sufficient to describe a binary float.
sdsf          := ceil( mantis_length * log10( 2 ) ) :

#-----

nearestBinaryFloat := proc( x :: numeric ) 
#Finds the binary float that is closest to x.

  local rd , xp , sg , mantis :

  #Special case
  if ( x = 0 ) then
    return 1 , emin - 1 , mantis_min :
  end if :

  #-----

  #When we apply Scale2, we need to compute the
  #integer part of the result exactly. So ...
  rd := ilog10( mantis_max ) + 5 : #We add a few extra digits.

  if ( Digits < rd ) then
    Digits := rd :
  end if :

  sg := sign( x ) :
  xp := ilog2( abs( x ) ) :

  if ( xp < emin ) then
    #return zero
    return 1 , emin - 1 , mantis_min :
  elif ( xp > emax ) then
    error "Result is not a machine number (overflow)" :
  end if :

  #Removing the exponent and dividing by 2 yields a number between
  #between 0 and 1. Multiplying by 2^mantis_length and rounding to
  #the nearest integer then yields the best possible mantissa.
  mantis := round( Scale2( abs( x ) , mantis_length - 1 - xp ) ) :

  #Note the trap here: if |x| * 2^(-1-xp) turns out to be too 
  #close to 1, then mantis will now be set to mantis_max + 1. So...
  if ( mantis > mantis_max ) then

    if ( xp = emax ) then
      error "Result is not a machine number (overflow)" :
    end if :

    mantis := mantis_min :
    xp     := xp + 1     :

  end if :

  #-----

  return sg , xp , mantis :

  #-----

end proc : 

#-----  

toFraction := proc( sg :: integer , xp :: integer , mantis :: integer ) 
#Converts a binary float to a fraction

  if ( xp < emin and mantis = mantis_min ) then
    return 0 :
  else

    #sg, xp and mantis are integers, so we need 
    #not worry about the value of Digits here.
    return sg * Scale2( mantis , 1 + xp - mantis_length ) :

  end if :

end proc :

#-----

stepBinaryFloat := proc( sg :: integer , xp :: integer , mantis :: integer , direction :: integer )
#If direction = 1, this function returns the smallest binary float that is larger than the number 
#represented by the ordered triple ( sg , mantis , xp). If direction = -1, it returns the largest 
#binary float that is smaller than ( sg , mantis , xp ).

  local sgd , x , m :

  #-----

  #Special case: input values represent zero.
  if ( xp < emin ) then
    return sign( direction ) , emin , mantis_min :
  end if :

  #-----

  sgd := sg * direction :

  if ( sgd = 1 ) then 
  #In this case we are increasing the number's magnitude

    if ( mantis = mantis_max ) then

      if ( xp = emax ) then 
      # x = emax + 1 is reserved for NaN, etc.
        error "Result is not a machine number (overflow)" :
      end if :

      m := mantis_min :
      x := xp + 1 :
 
    else

      m := mantis + 1 :
      x := xp :

    end if :

  else
  #In this case we are decreasing the number's magnitude

    if ( mantis = mantis_min ) then

      x := xp - 1 :

      if ( x < emin ) then 
        m := mantis_min : #Return zero
      else
        m := mantis_max :
      end if :

    else

      m := mantis - 1 :
      x := xp :

    end if :

  end if :

  #-----

  return sg , x , m :

  #-----

end proc :

#-----

transfer := proc( sg :: integer , xp :: integer , mantis :: integer )
#If x is the nbits bit real number represented by the triple ( sg , mantis , xp ),
#then the value returned by this function will be equal to that which would be 
#returned by the Fortran function transfer( x , 1_li ) where li is a kind type 
#parameter that specifies an nbits bit integer.

  local t , top_bit_val :

  #-----

  top_bit_val := 2^(nbits-1) :

  #The first term is the biased exponent, shifted up to occupy bits starting
  #with the second. The second term is the mantissa with the hidden bit removed. 
  t := Scale2( xp + bias , mantis_length - 1 ) + ( mantis - mantis_min ) :

  #Integers are stored in two's complement format,
  #so if we want to set the leading bit to 1...
  if ( sg = -1 ) then
    t := t - top_bit_val :
  end if :

  return t :

end proc :

#-----

LibraryTools[Save]( 'nbits' , 'mantis_length' , 'emax' , 'emin' , 'bias' , 
                    'mantis_min' , 'mantis_max' , 'sdsf' , "./" ) :
LibraryTools[Save]( 'nearestBinaryFloat' , 'toFraction' , 'stepBinaryFloat' , 'transfer' , "./" ) :

#-----

#End of file 'binary float.mpl'
