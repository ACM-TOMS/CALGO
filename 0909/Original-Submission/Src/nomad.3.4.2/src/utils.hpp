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
  \file   utils.hpp
  \brief  Utility functions (headers)
  \author Sebastien Le Digabel
  \date   2010-03-23
  \see    utils.cpp
*/
#ifndef __UTILS__
#define __UTILS__

#include <cctype>
#include <vector>
#include <set>
#include <list>
#include <iomanip>
#include <cmath>

// use of 'access' or '_access', and getpid() or _getpid() :
#ifdef _MSC_VER
#include <io.h>
#include <process.h>
#else
#include <unistd.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "defines.hpp"

namespace NOMAD {

  /// Convert a string into a NOMAD::bb_input_type.
  /**
     \param  s    The string               -- \b IN.
     \param  bbit The NOMAD::bb_input_type -- \b OUT.
     \return      A boolean equal to \c true if the conversion was possible.
  */
  bool string_to_bb_input_type ( const std::string & s , NOMAD::bb_input_type & bbit );

  /// Convert a string into a NOMAD::bb_output_type.
  /**
     \param s    The string                -- \b IN.
     \param bbot The NOMAD::bb_output_type -- \b OUT.
     \return     A boolean equal to \c true if the conversion was possible.
  */
  bool string_to_bb_output_type ( const std::string & s , NOMAD::bb_output_type & bbot );
  
  /// Convert a string into a NOMAD::hnorm_type.
  /**
     \param  s  The string            -- \b IN.
     \param  hn The NOMAD::hnorm_type -- \b OUT.
     \return    A boolean equal to \c true if the conversion was possible.
  */
  bool string_to_hnorm_type ( const std::string & s , NOMAD::hnorm_type & hn );

  /// Convert a string into a multi_formulation_type.
  /**
     \param  s   The string                        -- \b IN.
     \param  mft The NOMAD::multi_formulation_type -- \b OUT.
     \return     A boolean equal to \c true if the conversion was possible.
  */
  bool string_to_multi_formulation_type ( const std::string             & s   ,
					  NOMAD::multi_formulation_type & mft   );

  /// Convert a string with format "i-j" into two integers i and j.
  /**
     If s=="*" and if n is defined, then i=0 and j=*n-1.

     \param  s The string              -- \b IN.
     \param  i The first integer \c i  -- \b OUT.
     \param  j The second integer \c j -- \b OUT.
     \param  n Number of variables; use \c NULL if unknown
               -- \b IN -- \b optional (default = \c NULL).
     \param  check_order A boolean indicating if \c i and \c j are to be compared
             -- \b IN -- \b optional (default = \c true).
     \return A boolean equal to \c true if the conversion was possible.
  */
  bool string_to_index_range ( const std::string & s                  ,
			       int               & i                  ,
			       int               & j                  ,
			       int               * n           = NULL ,
			       bool                check_order = true   );

  /// Convert a string in {"YES","NO","Y","N"} to a boolean.
  /**
     \param s The string -- \b IN.
     \return  An integer equal to \c 0 for \c false, \c 1 for \c true,
              and \c -1 if the conversion failed.
  */
  int string_to_bool ( const std::string & s );
  
  
  /// Interpret a list of strings as a direction type.
  /**
     \param ls The list of strings -- \b IN.
     \param dt The NOMAD::direction_type -- \b OUT.
     \return   A boolean equal to \c true if the conversion was possible.
  */
  bool strings_to_direction_type ( const std::list<std::string> & ls ,
				   NOMAD::direction_type        & dt   );

  /// If a NOMAD::bb_output_type variable corresponds to a constraint.
  /**
     \param bbot The NOMAD::bb_output_type -- \b IN.
     \return     A boolean equal to \c true if \c bbot corresponds to a constraint.
  */
  bool bbot_is_constraint ( NOMAD::bb_output_type bbot );

  /// If a NOMAD::direction_type variable corresponds to a MADS direction.
  /**
     \param dt The NOMAD::direction_type -- \b IN.
     \return   A boolean equal to \c true if \c dt corresponds to a MADS direction.
  */
  bool dir_is_mads ( NOMAD::direction_type dt );

  /// If a NOMAD::direction_type variable corresponds to a GPS direction.
  /**
     \param dt The NOMAD::direction_type -- \b IN.
     \return   A boolean equal to \c true if \c dt corresponds to a GPS direction.
  */
  bool dir_is_gps ( NOMAD::direction_type dt );
  
  /// If a NOMAD::direction_type variable corresponds to a LT-MADS direction.
  /**
     \param dt The NOMAD::direction_type -- \b IN.
     \return   A boolean equal to \c true if \c dt corresponds to a LT-MADS direction.
  */
  bool dir_is_ltmads ( NOMAD::direction_type dt );

  /// If a NOMAD::direction_type variable corresponds to a random direction.
  /**
     \param dt The NOMAD::direction_type -- \b IN.
     \return   A boolean equal to \c true if \c dt corresponds to a random direction.
  */
  bool dir_is_random ( NOMAD::direction_type dt );


  /// If a NOMAD::direction_type variable corresponds to a Ortho-MADS direction.
  /**
     \param dt The NOMAD::direction_type -- \b IN.
     \return   A boolean equal to \c true if \c dt corresponds to a Ortho-MADS direction.
  */
  bool dir_is_orthomads ( NOMAD::direction_type dt );

  /// Check if a set of directions include Ortho-MADS direction.
  /**
     \param dir_types Set of direction types -- \b IN.
     \return A boolean equal to \c true if at
     least one direction in the set is
     of type Ortho-MADS.
  */
  bool dirs_are_orthomads ( const std::set<NOMAD::direction_type> & dir_types );

  /// Construct the n first prime numbers.
  /**
     \param n      The integer \c n-- \b IN.
     \param primes An integer array of size \c n for the prime numbers;
                   must be previously allocated -- \b OUT.
  */
  void construct_primes ( int n , int * primes );
  
  /// Decompose a string (sentence) into a list of strings (words).
  /**
     \param sentence The sentence -- \b IN.
     \param words    The words    -- \b OUT.
  */
  void get_words ( const std::string & sentence , std::list<std::string> & words );

  /// Check if a file exists and is executable.
  /**
     \param file_name A string corresponding to a file name -- \b IN.
     \return          A boolean equal to \c true if the file is executable.
  */
  bool check_exe_file  ( const std::string & file_name );

  /// Check if a file exists and is readable.
  /**
     \param file_name A string corresponding to a file name -- \b IN.
     \return          A boolean equal to \c true if the file exists and is readable.
  */
  bool check_read_file ( const std::string & file_name );

  /// Get the process id (pid); useful for unique random seeds.
  /**
     \return An integer corresponding to the pid.
  */
  int get_pid ( void );

  /// Called at the beginning of NOMAD.
  /**
     \param argc Number of command line arguments.
     \param argv Command line arguments.
  */
  void begin ( int argc , char ** argv );

  /// Called at the end of NOMAD.
  void end ( void );

  
  /// Transform an integer into a string.
  /**
     \param i The integer -- \b IN.
     \return  The string.
  */
  std::string itos ( int i );

  /// Put a string into upper cases.
  /**
     \param s The string -- \b IN/OUT.
  */
  void toupper ( std::string & s );

  /// Put a list of strings into upper cases.
  /**
     \param ls The list of strings -- \b IN/OUT.
  */
  void toupper  ( std::list<std::string> & ls );

  /// Convert a string into an integer.
  /**
     \param s The string  -- \b IN.
     \param i The integer -- \b OUT.
     \return  A boolean equal to \c true if the conversion was possible.
  */
  bool atoi ( const std::string & s , int & i );

  /// Convert a character into an integer.
  /**
     \param c The character -- \b IN.
     \param i The integer   -- \b OUT.
     \return  A boolean equal to \c true if the conversion was possible.
  */
  bool atoi ( char c , int & i );
 
  /// Search a list of string inside a string.
  /**
     \param  s  The string          -- \b IN.
     \param  ls The list of strings -- \b IN.
     \return    A boolean equal to \c true if one of the string of ls is in s.
  */
  bool string_find ( const std::string & s , const std::list<std::string> & ls );

  /// Search a string into another string.
  /**
     \param  s1 A string -- \b IN.
     \param  s2 A string -- \b IN.
     \return    A boolean equal to \c true if \c s2 is in \c s1.
  */
  bool string_find ( const std::string & s1 , const std::string & s2 );
}

#endif
