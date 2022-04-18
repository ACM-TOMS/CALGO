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
  \file   utils.cpp
  \brief  Utility functions
  \author Sebastien Le Digabel
  \date   2010-03-23
  \see    utils.hpp
*/
#include "utils.hpp"

/*---------------------------------------------------------------*/
/*               construct the first n prime numbers             */
/*---------------------------------------------------------------*/
void NOMAD::construct_primes ( int n , int * primes )
{
  bool   is_prime;
  double kk;
  int    i = 0 , k = 2 , j;
  while ( true ) {
    is_prime = true;
    kk = sqrt ( static_cast<double>(k) );
    for ( j = 2 ; j <= kk ; ++j )
      if ( k%j==0 ) {
	is_prime = false;
	break;
      }
    if ( is_prime ) {
      primes[i++] = k;
      if ( i==n )
	break;
    }
    ++k;
  }
}

/*--------------------------------------------*/
/*  decompose a string (sentence) into words  */
/*--------------------------------------------*/
void NOMAD::get_words ( const std::string & sentence , std::list<std::string> & words )
{
  std::string        s;
  std::istringstream in ( sentence );
   while ( true ) {
    in >> s;
    if ( in.fail() )
      break;
    words.push_back ( s );
  }
}

/*---------------------------------------------------------------*/
/*  get the pid (process id): used as random seed or unique tag  */
/*---------------------------------------------------------------*/
int NOMAD::get_pid ( void )
{
#ifdef _MSC_VER
  return _getpid();
#else
  return getpid();
#endif
}

/*---------------------------------------*/
/*    called at the beginning of NOMAD   */
/*---------------------------------------*/
void NOMAD::begin ( int argc , char ** argv )
{
#ifdef USE_MPI
  MPI_Init ( &argc , &argv );
#endif
}

/*---------------------------------------*/
/*      called at the end of NOMAD       */
/*---------------------------------------*/
void NOMAD::end ( void )
{
#ifdef USE_MPI
  MPI_Finalize();
#endif
}

/*-----------------------------------------------------------------*/
/*              check if a file exists and is executable           */
/*-----------------------------------------------------------------*/
bool NOMAD::check_exe_file ( const std::string & file_name )
{
#ifdef WINDOWS
  // don't check on Windows:
  return true;
#else
  return ( access ( file_name.c_str() , X_OK ) == 0 );
#endif
}

/*-----------------------------------------------------------------*/
/*              check if a file exists and is readable             */
/*-----------------------------------------------------------------*/
bool NOMAD::check_read_file ( const std::string & file_name )
{
#ifdef _MSC_VER
  return ( _access ( file_name.c_str() , 4 ) == 0 );
#else
  return ( access ( file_name.c_str() , R_OK ) == 0 );
#endif
}

/*-----------------------------------------------------------------*/
/*                         NOMAD::itos                             */
/*-----------------------------------------------------------------*/
std::string NOMAD::itos ( int i )
{
  std::ostringstream oss;
  oss << i;
  return oss.str();
}

/*-----------------------------------------------------------------*/
/*                         NOMAD::toupper - 1/2                    */
/*-----------------------------------------------------------------*/
void NOMAD::toupper ( std::string & s )
{
  size_t ns = s.size();
  for ( size_t i = 0 ; i < ns ; ++i )
    s[i] = std::toupper(s[i]);
}

/*-----------------------------------------------------------------*/
/*                         NOMAD::toupper - 2/2                    */
/*-----------------------------------------------------------------*/
void NOMAD::toupper ( std::list<std::string> & ls )
{
  std::list<std::string>::iterator       it;
  std::list<std::string>::const_iterator end = ls.end();
  for ( it = ls.begin() ; it != end ; ++it )
    NOMAD::toupper ( *it );
}

/*-----------------------------------------------------------------*/
/*                             NOMAD::atoi                         */
/*-----------------------------------------------------------------*/
bool NOMAD::atoi ( const std::string & s , int & i )
{
  i = -1;
  if ( s.empty() )
    return false;

  size_t n = s.size();

  if ( s[0] == '-' ) {
    if ( n > 1 && s[1] == '-' )
      return false;
    std::string ss = s;
    ss.erase(ss.begin());
    if ( NOMAD::atoi ( ss , i ) ) {
      i = -i;
      return true;
    }
    return false;
  }

  for ( size_t k = 0 ; k < n ; ++k )
    if ( !isdigit(s[k]) )
      return false;
  i = std::atoi(s.c_str());
  return true;
}

bool NOMAD::atoi ( char c , int & i ) {
  std::string s = "-";
  s[0] = c;
  return NOMAD::atoi(s,i);
}

/*-------------------------------------------------------------------*/
/*  if a NOMAD::bb_output_type variable corresponds to a constraint  */
/*-------------------------------------------------------------------*/
bool NOMAD::bbot_is_constraint ( NOMAD::bb_output_type bbot )
{
  return ( bbot == NOMAD::EB     ||
	   bbot == NOMAD::PB     ||
	   bbot == NOMAD::PEB_P  ||
	   bbot == NOMAD::PEB_E  ||
	   bbot == NOMAD::FILTER    );
}

/*-----------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a MADS direction  */
/*-----------------------------------------------------------------------*/
bool NOMAD::dir_is_mads ( NOMAD::direction_type dt )
{
  return ( dt == NOMAD::ORTHO_1  ||
	   dt == NOMAD::ORTHO_2  ||
	   dt == NOMAD::ORTHO_2N ||
	   dt == NOMAD::LT_1     ||
	   dt == NOMAD::LT_2     ||
	   dt == NOMAD::LT_2N    ||
	   dt == NOMAD::LT_NP1      );
}

/*----------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a GPS direction  */
/*----------------------------------------------------------------------*/
bool NOMAD::dir_is_gps ( NOMAD::direction_type dt )
{
  return ( dt == NOMAD::GPS_BINARY             ||
	   dt == NOMAD::GPS_2N_STATIC          ||
	   dt == NOMAD::GPS_2N_RAND            ||
	   dt == NOMAD::GPS_NP1_STATIC_UNIFORM ||
	   dt == NOMAD::GPS_NP1_STATIC         ||
	   dt == NOMAD::GPS_NP1_RAND_UNIFORM   ||
	   dt == NOMAD::GPS_NP1_RAND              );
}

/*--------------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a LT-MADS direction  */
/*--------------------------------------------------------------------------*/
bool NOMAD::dir_is_ltmads ( NOMAD::direction_type dt )
{
  return ( dt == NOMAD::LT_1     ||
	   dt == NOMAD::LT_2     ||
	   dt == NOMAD::LT_2N    ||
	   dt == NOMAD::LT_NP1      );
}

/*-------------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a random direction  */
/*-------------------------------------------------------------------------*/
bool NOMAD::dir_is_random ( NOMAD::direction_type dt )
{
  return ( dt == NOMAD::GPS_NP1_RAND         ||
	   dt == NOMAD::GPS_NP1_RAND_UNIFORM ||
	   dt == NOMAD::GPS_2N_RAND             );
}

/*-----------------------------------------------------------------------------*/
/*  if a NOMAD::direction_type variable corresponds to a Ortho-MADS direction  */
/*-----------------------------------------------------------------------------*/
bool NOMAD::dir_is_orthomads ( NOMAD::direction_type dt )
{
  return ( dt == NOMAD::ORTHO_1  ||
	   dt == NOMAD::ORTHO_2  ||
	   dt == NOMAD::ORTHO_2N    );
}

/*---------------------------------------------------------------------*/
/*  check if a set of directions include Ortho-MADS direction          */
/*  (true if at least one direction in the set is of type Ortho-MADS)  */
/*---------------------------------------------------------------------*/
bool NOMAD::dirs_are_orthomads ( const std::set<NOMAD::direction_type> & dir_types )
{
  std::set<NOMAD::direction_type>::const_iterator it , end = dir_types.end();
  for ( it = dir_types.begin() ; it != end ; ++it )
    if ( NOMAD::dir_is_orthomads (*it) )
      return true;
  return false;
}

/*---------------------------------------------------*/
/*  returns true if one of the string of ls is in s  */
/*---------------------------------------------------*/
bool NOMAD::string_find ( const std::string & s , const std::list<std::string> & ls )
{
  std::list<std::string>::const_iterator it , end = ls.end();
  for ( it = ls.begin() ; it != end ; ++it )
    if ( NOMAD::string_find ( s , *it ) )
      return true;
  return false;
}

/*---------------------------------------------------*/
/*        search a string into another string        */
/*---------------------------------------------------*/
bool NOMAD::string_find ( const std::string & s1 , const std::string & s2 )
{
  return ( s1.find(s2) < s1.size() );
}

/*-----------------------------------------------------------------*/
/*         interpret a list of strings as a direction type         */
/*-----------------------------------------------------------------*/
bool NOMAD::strings_to_direction_type ( const std::list<std::string> & ls ,
					NOMAD::direction_type        & dt   )
{
  
  dt = NOMAD::UNDEFINED_DIRECTION;

  if ( ls.empty() || ls.size() > 4 )
    return false;

  std::list<std::string>::const_iterator it = ls.begin() , end = ls.end();
  std::string                            s  = *it;
  NOMAD::toupper ( s );

  // no direction:
  if ( s == "NONE" ) {
    dt = NOMAD::NO_DIRECTION;
    return true;
  }

  // Ortho-MADS with 1, 2 or 2n directions:
  if ( s == "ORTHO" ) {
    ++it;
    if ( it == end ) {
      dt = NOMAD::ORTHO_2N;
      return true;
    }
    if ( *it == "1" ) {
      dt = NOMAD::ORTHO_1;
      return true;
    }
    if ( *it == "2" ) {
      dt = NOMAD::ORTHO_2;
      return true;
    }
    s = *it;
    NOMAD::toupper ( s );
    if ( s == "2N" ) {
      dt = NOMAD::ORTHO_2N;
      return true;
    }
    return false;
  }

  // LT-MADS with 1, 2 or 2n directions:
  if ( s == "LT" ) {
    ++it;
    if ( it == end ) {
      dt = NOMAD::LT_2N;
      return true;
    }
    if ( *it == "1" ) {
      dt = NOMAD::LT_1;
      return true;
    }
    if ( *it == "2" ) {
      dt = NOMAD::LT_2;
      return true;
    }
    s = *it;
    NOMAD::toupper ( s );
    if ( s == "N+1" ) {
      dt = NOMAD::LT_NP1;
      return true;
    }
    if ( s == "2N" ) {
      dt = NOMAD::LT_2N;
      return true;
    }
    return false;
  }

  // GPS:
  if ( s == "GPS" ) {
    ++it;
    if ( it == end ) {
      dt = NOMAD::GPS_2N_STATIC;
      return true;
    }
    s = *it;
    NOMAD::toupper ( s );

    // GPS for binary variables:
    if ( s == "BINARY" || s == "BIN" ) {
      dt = NOMAD::GPS_BINARY;
      return true;
    }

    // GPS, n+1 directions:
    if ( s == "N+1" ) {
      ++it;
      if ( it == end ) {
	dt = NOMAD::GPS_NP1_STATIC;
	return true;
      }
      s = *it;
      NOMAD::toupper ( s );

      // GPS, n+1, static:
      if ( s == "STATIC" ) {
	++it;
	if ( it == end ) {
	  dt = NOMAD::GPS_NP1_STATIC;
	  return true;
	}
	s = *it;
	NOMAD::toupper ( s );
	if ( s == "UNIFORM" ) {
	  dt = NOMAD::GPS_NP1_STATIC_UNIFORM;
	  return true;
	}
	return false;
      }

      // GPS, n+1, random:
      if ( s == "RAND" || s == "RANDOM" ) {
	++it;
	if ( it == end ) {
	  dt = NOMAD::GPS_NP1_RAND;
	  return true;
	}
	s = *it;
	NOMAD::toupper ( s );
	if ( s == "UNIFORM" ) {
	  dt = NOMAD::GPS_NP1_RAND_UNIFORM;
	  return true;
	}
	return false;
      }
      return false;
    }
    
    // 2n directions:
    if ( s == "2N" ) {
      ++it;
      if ( it == end ) {
	dt = NOMAD::GPS_2N_STATIC;
	return true;
      }
      s = *it;
      NOMAD::toupper ( s );
      if ( s == "STATIC" ) {
	dt = NOMAD::GPS_2N_STATIC;
	return true;
      }
      if ( s == "RAND" || s == "RANDOM" ) {
	dt = NOMAD::GPS_2N_RAND;
	return true;
      }
      return false;
    }
    return false;
  }
  return false;
}

/*---------------------------------------*/
/*   convert a string into a hnorm_type  */
/*---------------------------------------*/
bool NOMAD::string_to_hnorm_type ( const std::string & s , NOMAD::hnorm_type & hn )
{
  std::string ss = s;
  NOMAD::toupper(ss);
  if ( ss == "L1" ) {
    hn = NOMAD::L1;
    return true;
  }
  if ( ss == "L2" ) {
    hn = NOMAD::L2;
    return true;
  }
  if ( ss == "LINF" ) {
    hn = NOMAD::LINF;
    return true;
  }
  return false;
}

/*--------------------------------------------------*/
/*  convert a string into a multi_formulation_type  */
/*--------------------------------------------------*/
bool NOMAD::string_to_multi_formulation_type ( const std::string             & s   ,
					       NOMAD::multi_formulation_type & mft   )
{
  std::string ss = s;
  NOMAD::toupper(ss);
  if ( ss == "NORMALIZED" ) {
    mft = NOMAD::NORMALIZED;
    return true;
  }
  if ( ss == "PRODUCT" ) {
    mft = NOMAD::PRODUCT;
    return true;
  }
  if ( ss == "DIST_L1" ) {
    mft = NOMAD::DIST_L1;
    return true;
  }
  if ( ss == "DIST_L2" ) {
    mft = NOMAD::DIST_L2;
    return true;
  }
  if ( ss == "DIST_LINF" ) {
    mft = NOMAD::DIST_LINF;
    return true;
  }
  return false;
}

/*-------------------------------------------*/
/*   convert a string into a bb_output_type  */
/*-------------------------------------------*/
bool NOMAD::string_to_bb_output_type ( const std::string     & s    ,
				       NOMAD::bb_output_type & bbot   )
{
  std::string ss = s;
  NOMAD::toupper(ss);
  
  if ( ss == "OBJ" ) {
    bbot = NOMAD::OBJ;
    return true;
  }
  if ( ss == "EB" ) {
    bbot = NOMAD::EB;
    return true;
  }
  if ( ss == "PB" || ss == "CSTR" ) {
    bbot = NOMAD::PB;
    return true;
  }
  if ( ss == "PEB" ) {
    bbot = NOMAD::PEB_P;
    return true;
  }
  if ( ss == "F" ) {
    bbot = NOMAD::FILTER;
    return true;
  }
  if ( ss == "STAT_AVG" ) {
    bbot = NOMAD::STAT_AVG;
    return true;
  }
  if ( ss == "STAT_SUM" ) {
    bbot = NOMAD::STAT_SUM;
    return true;
  }
  if ( ss == "CNT_EVAL" ) {
    bbot = NOMAD::CNT_EVAL;
    return true;
  }
  if ( ss == "NOTHING" || ss == "-" ) {
    bbot = NOMAD::UNDEFINED_BBO;
    return true;
  }
  return false;
}

/*-----------------------------------------------------------------*/
/*                convert a string into a bb_input_type            */
/*-----------------------------------------------------------------*/
bool NOMAD::string_to_bb_input_type ( const std::string    & s    ,
				      NOMAD::bb_input_type & bbit   )
{
  std::string ss = s;
  NOMAD::toupper ( ss );
  if ( ss=="R" || ss=="REAL" ) {
    bbit = NOMAD::CONTINUOUS;
    return true;
  }
  if ( ss=="C" || ss=="CAT" ) {
    bbit = NOMAD::CATEGORICAL;
    return true;
  }
  if ( ss=="B" || ss=="BIN" ) {
    bbit = NOMAD::BINARY;
    return true;
  }
  if ( ss=="I" || ss=="INT" ) {
    bbit = NOMAD::INTEGER;
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------*/
/*         convert a string in {"YES","NO","Y","N"} to a bool           */
/*         value of return: -1: error                                   */
/*                           0: bool=false                              */
/*                           1: bool=true                               */
/*----------------------------------------------------------------------*/
int NOMAD::string_to_bool ( const std::string & ss )
{
  std::string s = ss;
  NOMAD::toupper ( s );
  if ( s=="Y" || s=="YES" || s=="1" )
    return 1;
  if ( s=="N" || s=="NO" || s=="0" )
    return 0;
  return -1;
}

/*---------------------------------------------------*/
/*  convert a string 'i-j' to the integers i and j   */
/*---------------------------------------------------*/
bool NOMAD::string_to_index_range ( const std::string & s           ,
				    int               & i           ,
				    int               & j           ,
				    int               * n           ,
				    bool                check_order   )
{
  if ( s.empty() )
    return false;
  
  if ( s == "*" ) {
    if ( !n )
      return false;
    i = 0;
    j = *n-1;
    return true;
  }
  
  if ( s[0] == '-' ) {
      
    size_t ns = s.size();
    if ( ns > 1 && s[1] == '-' )
      return false;
    
    std::string ss = s;
    ss.erase ( ss.begin() );
      
    if ( NOMAD::string_to_index_range ( ss , i , j , n , false ) ) {
      i = -i;
      return true;
    }
    return false;
  }

  std::istringstream in (s);
  std::string        s1;

  getline ( in , s1 , '-' );

  if (in.fail())
    return false;

  size_t k , n1 = s1.size();

  if ( n1 >= s.size() - 1 ) {
    for ( k = 0 ; k < n1 ; ++k )
      if (!isdigit(s1[k]))
	return false;
    if ( !NOMAD::atoi ( s1 , i ) )
      return false;
    if ( n1 == s.size() ) {
      j = i;
      return true;
    }
    if (n) {
      j = *n-1;
      return true;
    }
    return false;
  }

  std::string s2;
  getline (in, s2);

  if (in.fail())
    return false;

  size_t n2 = s2.size();
  for ( k = 0 ; k < n2 ; ++k )
    if ( !isdigit(s2[k]) )
      return false;

  if ( !NOMAD::atoi ( s1, i ) || !NOMAD::atoi ( s2 , j ) )
    return false;

  return !check_order || i <= j;
}
