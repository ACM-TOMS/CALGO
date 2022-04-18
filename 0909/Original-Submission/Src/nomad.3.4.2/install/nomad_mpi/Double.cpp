/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonsmooth Optimization by Mesh Adaptive Direct search - version 3.4        */
/*                                                                                     */
/*  Copyright (C) 2001-2010  Mark Abramson        - the Boeing Company, Seattle        */
/*                           Charles Audet        - Ecole Polytechnique, Montreal      */
/*                           Gilles Couture       - Ecole Polytechnique, Montreal      */
/*                           John Dennis          - Rice University, Houston           */
/*                           Sebastien Le Digabel - Ecole Polytechnique, Montreal      */
/*                                                                                     */
/*  funded in part by AFOSR and Exxon Mobil                                            */
/*                                                                                     */
/*  Author: Sebastien Le Digabel                                                       */
/*                                                                                     */
/*  Contact information:                                                               */
/*    Ecole Polytechnique de Montreal - GERAD                                          */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada                  */
/*    e-mail: nomad@gerad.ca                                                           */
/*    phone : 1-514-340-6053 #6928                                                     */
/*    fax   : 1-514-340-5665                                                           */
/*                                                                                     */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad               */
/*-------------------------------------------------------------------------------------*/
/**
  \file   Double.cpp
  \brief  Custom class for double-precision reals (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-02
  \see    Double.hpp
*/
#include "Double.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
double      NOMAD::Double::_epsilon         = NOMAD::DEFAULT_EPSILON;
std::string NOMAD::Double::_inf_str         = NOMAD::DEFAULT_INF_STR;
std::string NOMAD::Double::_undef_str       = NOMAD::DEFAULT_UNDEF_STR;
#ifdef DEBUG
int         NOMAD::Double::_cardinality     = 0;
int         NOMAD::Double::_max_cardinality = 0;
#endif

/*-----------------------------------------------*/
/*                  constructor 1                */
/*-----------------------------------------------*/
NOMAD::Double::Double ( void )
  : _value   ( 0.0   ) ,
    _defined ( false )
{
#ifdef DEBUG
  ++NOMAD::Double::_cardinality;
  if ( NOMAD::Double::_cardinality > NOMAD::Double::_max_cardinality )
    ++NOMAD::Double::_max_cardinality;
#endif
}

/*-----------------------------------------------*/
/*                  constructor 2                */
/*-----------------------------------------------*/
NOMAD::Double::Double ( double v )
  : _value   ( v    ) ,
    _defined ( true )
{
#ifdef DEBUG
  ++NOMAD::Double::_cardinality;
  if (NOMAD::Double::_cardinality > NOMAD::Double::_max_cardinality)
    ++NOMAD::Double::_max_cardinality;
#endif
}

/*-----------------------------------------------*/
/*                  constructor 3                */
/*-----------------------------------------------*/
NOMAD::Double::Double ( const NOMAD::Double & d )
  : _value   ( d._value   ) ,
    _defined ( d._defined )
{
#ifdef DEBUG
  ++NOMAD::Double::_cardinality;
  if (NOMAD::Double::_cardinality > NOMAD::Double::_max_cardinality)
    ++NOMAD::Double::_max_cardinality;
#endif
}

/*-----------------------------------------------*/
/*                    destructor                 */
/*-----------------------------------------------*/
NOMAD::Double::~Double ( void ) 
{
#ifdef DEBUG
  _value   = 0.0;
  _defined = false;
  --NOMAD::Double::_cardinality;
#endif
}

/*-----------------------------------------------*/
/*               set epsilon (static)            */
/*-----------------------------------------------*/
void NOMAD::Double::set_epsilon ( double eps ) 
{
  if ( eps <= 0.0 )
    throw NOMAD::Exception ( "Double.cpp" , __LINE__ ,
			     "NOMAD::Double::set_epsilon(): invalid epsilon" );
  NOMAD::Double::_epsilon = eps;
}

/*-----------------------------------------------*/
/*                  get the value                */
/*-----------------------------------------------*/
double NOMAD::Double::value ( void ) const
{
  if ( !_defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double::value(): value not defined" );
  return _value;
}

/*------------------------------------------*/
/*                    input                 */
/*------------------------------------------*/
std::istream & NOMAD::operator >> ( std::istream & in , NOMAD::Double & d )
{
  std::string s;
  in >> s;

  if ( !in.fail() && !d.atof (s) )    
    in.setstate ( std::ios::failbit );

  return in;
}

/*-----------------------------------------------*/
/*      atof: value determined by a string       */
/*-----------------------------------------------*/
bool NOMAD::Double::atof ( const std::string & ss )
{

  std::string s = ss;
  NOMAD::toupper(s);

  if ( s == "-" || ss == NOMAD::Double::_undef_str ) {
    _value   = 0.0;
    _defined = false;
    return true;
  }

  if ( s == "INF" ||  s == "+INF" ||
       ss == NOMAD::Double::_inf_str ||
       ss == ("+" + NOMAD::Double::_inf_str) ) {
    _value   = NOMAD::INF;
    _defined = true;
    return true;
  }

  if ( s == "-INF" || ss == ("-" + NOMAD::Double::_inf_str) ) {
    _value   = -NOMAD::INF;
    _defined = true;
    return true;
  }

  if ( s.empty() || (s.size() == 1 && !isdigit(s[0])) )
    return false;

  if ( !isdigit(s[0]) && s[0] != '+' && s[0] != '-' && s[0] != '.' )
    return false;

  size_t n = s.size();
  for ( size_t k = 1 ; k < n ; ++k )
    if ( !isdigit(s[k]) && s[k] != '.' ) {
      if ( s[k] == 'E' ) {
 	if ( s.size() == k+1 )
 	  return false;
 	++k;
 	if ( !isdigit(s[k]) && s[k] != '+' && s[k] != '-' )
 	  return false;
      }
      else
 	return false;
    }

  *this = std::atof ( s.c_str() );
  return true;
}

/*-------------------------------------------------*/
/*  atof from a string that can begin with 'r' to  */
/*  indicate a proportion (relative value)         */
/*-------------------------------------------------*/
bool NOMAD::Double::relative_atof ( const std::string & s , bool & relative )
{
  if ( std::toupper(s[0]) == 'R' ) {
    relative  = true;
    std::string ss = s;
    ss.erase(ss.begin());
    if ( !atof(ss) )
      return false;
    return ( *this >= 0.0 );
  }
  relative = false;
  return atof(s);
}

/*-----------------------------------------------*/
/*            is the value an integer?           */
/*-----------------------------------------------*/
bool NOMAD::Double::is_integer ( void ) const
{
  if ( !_defined )
    return false;
  return ( NOMAD::Double(floor(_value))) == ( NOMAD::Double(ceil(_value)) );
}

/*-----------------------------------------------*/
/*             is the value binary ?             */
/*-----------------------------------------------*/
bool NOMAD::Double::is_binary ( void ) const
{
  if ( !_defined )
    return false;
  return ( NOMAD::Double(_value) == 0.0 || NOMAD::Double(_value) == 1.0 );
}

/*-------------------------------------*/
/*               d = d1/d2             */
/*-------------------------------------*/
const NOMAD::Double NOMAD::operator / ( const NOMAD::Double & d1 ,
					const NOMAD::Double & d2   )
{
  if ( d2.value() == 0.0 )
    throw NOMAD::Double::Invalid_Value ( "Double.cpp" , __LINE__ ,
					 "NOMAD::Double: d1 / d2: division by zero" );
  return NOMAD::Double ( d1.value() / d2.value() );
}

/*-------------------------------------*/
/*                d1 += d2             */
/*-------------------------------------*/
const NOMAD::Double & NOMAD::Double::operator += ( const NOMAD::Double & d2 )
{
  if ( !_defined || !d2._defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double: d1 += d2: d1 or d2 not defined" );
  _value += d2._value;
  return *this;
}

/*-------------------------------------*/
/*               d1 -= d2              */
/*-------------------------------------*/
const NOMAD::Double & NOMAD::Double::operator -= ( const NOMAD::Double & d2 )
{
  if ( !_defined || !d2._defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double: d1 -= d2: d1 or d2 not defined" );
  _value -= d2._value;
  return *this;
}

/*-------------------------------------*/
/*                d1 *= d2             */
/*-------------------------------------*/
const NOMAD::Double & NOMAD::Double::operator *= ( const NOMAD::Double & d2 )
{
  if ( !_defined || !d2._defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double: d1 *= d2: d1 or d2 not defined" );
  _value *= d2._value;
  return *this;
}

/*-------------------------------------*/
/*               d1 /= d2              */
/*-------------------------------------*/
const NOMAD::Double & NOMAD::Double::operator /= ( const NOMAD::Double & d2 )
{
  if ( !_defined || !d2._defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double: d1 /= d2: d1 or d2 not defined" );
  if ( d2._value == 0.0 )
    throw Invalid_Value ( "Double.cpp" , __LINE__ ,
			  "NOMAD::Double: d1 /= d2: division by zero" );
  _value /= d2._value;
  return *this;
}

/*-------------------------------------*/
/*                  ++d                */
/*-------------------------------------*/
NOMAD::Double & NOMAD::Double::operator++ ( void )
{
  if ( !_defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ , "NOMAD::Double: ++d: d not defined" );
  _value += 1;
  return *this;
}

/*-------------------------------------*/
/*                  d++                */
/*-------------------------------------*/
NOMAD::Double NOMAD::Double::operator++ ( int n )
{
  if ( !_defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ , "NOMAD::Double: d++: d not defined" );
  NOMAD::Double tmp = *this;
  if( n <= 0 )
    n = 1;
  _value += n;
  return tmp;
}

/*-------------------------------------*/
/*                --d                  */
/*-------------------------------------*/
NOMAD::Double & NOMAD::Double::operator-- ( void )
{
  if ( !_defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ , "NOMAD::Double: --d: d not defined" );
  _value -= 1;
  return *this;
}

/*-------------------------------------*/
/*                  d--                */
/*-------------------------------------*/
NOMAD::Double NOMAD::Double::operator-- ( int n )
{
  if ( !_defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double: d--: d not defined" );
  NOMAD::Double tmp = *this;
  if ( n <= 0 )
    n = 1;
  _value -= n;
  return tmp;
}

/*-------------------------------------*/
/*              operators =            */
/*-------------------------------------*/
NOMAD::Double & NOMAD::Double::operator = ( const NOMAD::Double & d )
{
  _value   = d._value;
  _defined = d._defined;
  return *this;
}

NOMAD::Double & NOMAD::Double::operator = ( double r )
{
  _value   = r;
  _defined = true;
  return *this;
}

/*------------------------------------------*/
/*                  display                 */
/*------------------------------------------*/
void NOMAD::Double::display ( const NOMAD::Display & out ) const
{
  if ( _defined ) {
    if ( _value == NOMAD::INF )
      out << NOMAD::Double::_inf_str;
    else if ( _value == -NOMAD::INF )
      out << "-" << NOMAD::Double::_inf_str;
    else if ( floor(_value) == ceil(_value) && fabs(_value) < INT_MAX-1 )
      out << static_cast<int>(_value);
    else
      out << _value;
  }
  else
    out << NOMAD::Double::_undef_str;
}

/*------------------------------------------*/
/*              display with format         */
/*------------------------------------------*/
void NOMAD::Double::display ( const NOMAD::Display & out    ,
			      const std::string    & format   ) const
{

  // interpret the format:
  // ---------------------
 
  // %f      w=-1 prec=-1 c='f'
  // %4.5f   w= 4 prec= 5 c='f'
  // %4f     w= 4 prec= 1 c='f'
  // %.5f    w=-1 prec= 5 c='f'
  // %.f     w=-1 prec= 0 c='f'

  // c may be in 'e', 'E', 'f', 'g', 'G', 'd', or 'i'

  // e Scientific notation (mantise/exponent) using e character 3.9265e+2
  // E Scientific notation (mantise/exponent) using E character 3.9265E+2
  // f Decimal floating point                                   392.65
  // g Use the shorter of %e or %f                              392.65
  // G Use the shorter of %E or %f                              392.65
  // d or i Integer rounded value                               393

  std::string format2 = format;

  int  w    = -1;
  int  prec = -1;
  char c    =  0;

  if ( !format2.empty() && format2[0]=='%' ) {

    size_t n = format2.size();

    c = format2[n-1];
    
    if ( c!='e' && c!='E' && c!='f' && c!='g' && c!='G' && c!='d' && c!='i' ) {
      c = ( floor(_value) == ceil(_value) && fabs(_value) < INT_MAX-1 ) ?
	'd' : 'f';
      format2.push_back(c);
      ++n;
    }

    if ( n > 2 ) {

      std::string sw , sprec;

      size_t k = format2.find(".");
      
      if ( k > 0 && k < n-1 ) {
	if ( n==3 ) {
	  sprec = "0";
	}
	else {
	  if ( k > 1 )
	    sw = format2.substr ( 1 , k-1 );
	  sprec = format2.substr ( k+1 , n-k-2 );
	}
      }
      else {
	sw = format2.substr ( 1 , n-2 );
      }
      
      if ( !NOMAD::atoi ( sw , w ) )
	w = -1;

      if ( !NOMAD::atoi ( sprec , prec ) )
	prec = -1;
    }

    if ( c=='d' || c=='i' )
      prec = 0;
  }

  // display the value:
  out << std::setw(w);
  if ( _defined ) {
    if ( _value == NOMAD::INF )
      out << NOMAD::Double::_inf_str;
    else if ( c=='d' || c=='i' ||
	      ( format2.empty() &&
		floor(_value) == ceil(_value) && fabs(_value) < INT_MAX-1 ) )
      out << round();
    else {
      
      int                old_prec  = out.precision();
      std::_Ios_Fmtflags old_flags = out.flags();

      if ( prec >= 0 )
	out.precision ( prec );

      if ( c == 'f' )
	out.setf ( std::ios::fixed );

      else if ( c == 'e' ) {
	out.unsetf ( std::ios::fixed );
	out.setf   ( std::ios::scientific );
      }

      else if ( c == 'E' ) {
	out.unsetf ( std::ios::fixed );
	out.setf   ( std::ios::scientific | std::ios::uppercase );
      }

      else if ( c == 'g' ) {
	out.unsetf ( std::ios::fixed );
	out.setf   ( std::ios::floatfield );
      }

      else if ( c == 'G' ) {
	out.unsetf ( std::ios::fixed );
	out.setf   ( std::ios::floatfield | std::ios::uppercase );
      }

      out << _value;

      out.precision ( old_prec  );
      out.flags     ( old_flags );
    }
  }
  else
    out << NOMAD::Double::_undef_str;
}

/*------------------------------------------*/
/*                    round                 */
/*------------------------------------------*/
int NOMAD::Double::round ( void ) const
{
  if ( !_defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double::round(): value not defined" );
  return static_cast<int>(_value < 0.0 ? -floor(fabs(_value) + .5) : floor(_value + .5));
}

/*------------------------------------------*/
/*                    abs                   */
/*------------------------------------------*/
const NOMAD::Double NOMAD::Double::abs ( void ) const
{
  if ( !_defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double::abs(): value not defined" );
  return fabs ( _value );
}

/*------------------------------------------*/
/*                  square                  */
/*------------------------------------------*/
const NOMAD::Double NOMAD::Double::pow2 ( void ) const
{
  if ( !_defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double::pow2(): value not defined" );
  return pow ( _value , 2 );
}

/*------------------------------------------*/
/*                square root               */
/*------------------------------------------*/
const NOMAD::Double NOMAD::Double::sqrt ( void ) const
{
  if ( !_defined )
    throw Not_Defined ( "Double.cpp" , __LINE__ ,
			"NOMAD::Double::sqrt(): value not defined" );
  if ( *this < 0.0 )
    throw NOMAD::Double::Invalid_Value ( "Double.cpp" , __LINE__ ,
		        "NOMAD::Double::sqrt(x): x < 0" );

  return std::sqrt ( _value );
}

/*---------------------------------------------------------------------*/
/*  the same as operator < but with consideration of undefined values  */
/*---------------------------------------------------------------------*/
bool NOMAD::Double::comp_with_undef ( const NOMAD::Double & d ) const
{
  if ( this == &d )
    return false;

  bool d1d = is_defined();
  bool d2d = d.is_defined();

  if ( !d1d && !d2d )
    return false;

  if ( !d1d )
    return true;

  if ( !d2d )
    return false;

  return ( *this < d );
}

/*------------------------------------*/
/*  projection to mesh of size delta  */
/*  ( *this = ref + k * delta )       */
/*------------------------------------*/
void NOMAD::Double::project_to_mesh ( const NOMAD::Double & ref   ,
				      const NOMAD::Double & delta ,
				      const NOMAD::Double & lb    ,
				      const NOMAD::Double & ub      )
{
  if ( !_defined )
    return;

  NOMAD::Double v0 = ( ref._defined ) ? ref : v0;
  
  if ( delta._defined && delta != 0.0 ) {

    *this = v0 + ( (*this-v0) / delta).round() * delta;

    if ( ub._defined && *this > ub )
      *this = ub;

    if ( lb._defined && *this < lb )
      *this = lb;
  }
}
