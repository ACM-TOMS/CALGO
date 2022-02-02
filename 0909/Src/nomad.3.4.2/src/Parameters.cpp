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
  \file   Parameters.cpp
  \brief  NOMAD parameters (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-20
  \see    Parameters.hpp
*/
#include "Parameters.hpp"
#include "Slave.hpp"

/*----------------------------------------*/
/*                destructor              */
/*----------------------------------------*/
NOMAD::Parameters::~Parameters ( void )
{
  delete _std_signature;
  delete_x0s();
  reset_variable_groups();
}

/*----------------------------------------*/
/*         initializations (private)      */
/*----------------------------------------*/
void NOMAD::Parameters::init ( void )
{
  // miscellaneous and algorithm parameters:
  _to_be_checked      = true;
  _seed               = NOMAD::get_pid();
  _max_eval           = -1;
  _max_sim_bb_eval    = -1;
  _max_bb_eval        = -1;
  _max_bbe_decided    = false;
  _max_time           = -1;
  _max_iterations     = -1;
  _max_cache_memory   = -1.0;
  _cache_save_period  = 25;
  _snap_to_bounds     = true;
  _stop_if_feasible   = false;
  _user_calls_enabled = true;
  _asynchronous       = true;
  _stat_sum_target.clear();
  _L_curve_target.clear();
  _cache_file.clear();
  _problem_dir.clear();
  _tmp_dir.clear();

  // F_TARGET:
  reset_f_target();

  // output files:
  _add_seed_to_file_names = true;
  _solution_file.clear();
  _history_file.clear();

  // NOMAD::Double static members:
  set_EPSILON   ( NOMAD::DEFAULT_EPSILON   );
  set_UNDEF_STR ( NOMAD::DEFAULT_UNDEF_STR );
  set_INF_STR   ( NOMAD::DEFAULT_INF_STR   );

  // Mesh:
  _mesh_update_basis.clear();
  _mesh_coarsening_exponent = 1;
  _mesh_refining_exponent   =-1;
  _initial_mesh_index       = 0;
  _max_mesh_index           = NOMAD::UNDEFINED_L;
  _min_poll_size_defined    = false;
  _initial_mesh_size.clear();
  _min_mesh_size.clear();
  _min_poll_size.clear();
  
  // Directions:
  reset_directions ( -1 );

  // X0:
  reset_X0();

  // signature:
  delete _std_signature;
  _std_signature = NULL;  

  // variables:
  _dimension             = -1;
  _nb_free_variables     = -1;
  _extended_poll_trigger = 0.1;
  _relative_ept          = true;
  _extended_poll_enabled = true;
  _bb_input_type.clear();
  _bb_input_include_tag  = false;
  _bb_input_include_seed = false;
  _bb_redirection        = true;
  reset_bounds();
  reset_scaling();
  reset_fixed_variables();
  reset_periodic_variables();
  reset_variable_groups();
  
  // Barrier:
  _rho                    = 0.1;
  _h_min                  = 0.0;
  _h_max_0                = NOMAD::INF;
  _h_norm                 = NOMAD::L2;
  _has_constraints        = false;
  _has_EB_constraints     = false;
  _has_filter_constraints = false;
  _barrier_type           = NOMAD::EB;

  // outputs:
  _index_obj.clear();
  _bb_exe.clear();
  _bb_output_type.clear();
  _index_cnt_eval = -1;
  _index_stat_sum = -1;
  _index_stat_avg = -1;

  // surrogate:
  _sgte_exe.clear();
  _sgte_cache_file.clear();
  _sgte_eval_sort = true;
  _has_sgte       = false;
  _opt_only_sgte  = false;
  _sgte_cost      = -1;
  _sgte_max_eval  = -1;

  // MULTI-MADS:
  _multi_nb_mads_runs    = -1;
  _multi_overall_bb_eval = -1;
  _multi_use_delta_crit  = false;
  _multi_formulation     = NOMAD::UNDEFINED_FORMULATION;
  _multi_f_bounds.reset();

  // searches:
  _VNS_trigger.clear();
  _speculative_search         = true;
  _VNS_search                 = false;
  _LH_search_p0               = -1;
  _LH_search_pi               = -1;
  _opportunistic_LH           = true;
  _opp_LH_is_defined          = false;
  _cache_search               = false;
  _opportunistic_cache_search = false;
  _opp_CS_is_defined          = false;

  // opportunistic strategy:
  _opportunistic_eval           = true;
  _opportunistic_min_nb_success = -1;
  _opportunistic_min_eval       = -1;
  _opportunistic_lucky_eval     = false;
  _opportunistic_min_f_imprvmt.clear();

  // DISPLAY_DEGREE's set to max value for DEBUG:
#ifdef DEBUG
  _out.set_degrees ( NOMAD::FULL_DISPLAY );
  set_POINT_DISPLAY_LIMIT ( -1 );
#else
  _out.set_degrees ( NOMAD::NORMAL_DISPLAY );
  set_POINT_DISPLAY_LIMIT ( NOMAD::DEFAULT_POINT_DISPLAY_LIMIT );
#endif

  _out.set_open_brace   ( "{" );
  _out.set_closed_brace ( "}" );

  _display_stats.clear();
  reset_stats_file();
}

/*------------------------------------------------------------------*/
/*  delete the list of x0 points: std::vector<NOMAD::Point *> _x0s  */
/*  (private)                                                       */
/*------------------------------------------------------------------*/
void NOMAD::Parameters::delete_x0s ( void )
{
  size_t x0n = _x0s.size();
  for ( size_t i = 0 ; i < x0n ; ++i )
    delete _x0s[i];
  _x0s.clear();
}

/*----------------------------------------*/
/*           reset parameter X0           */
/*----------------------------------------*/
void NOMAD::Parameters::reset_X0 ( void )
{
  _to_be_checked = true;
  delete_x0s();
  _x0_cache_file.clear();
}

/*----------------------------------------*/
/*           reset the directions         */
/*----------------------------------------*/
void NOMAD::Parameters::reset_directions ( int halton_seed )
{
  _to_be_checked = true;
  _direction_types.clear();
  _sec_poll_dir_types.clear();
  _halton_seed = halton_seed;
}

/*--------------------------------------------*/
/*       check the display and file stats     */
/*--------------------------------------------*/
bool NOMAD::Parameters::check_display_stats
( const std::list<std::string> & stats ) const
{
  int var_index;
  std::list<std::string>::const_iterator it , end = stats.end();
  for ( it = stats.begin() ; it != end ; ++it ) {
    if ( !it->empty() &&
	 NOMAD::Display::get_display_stats_type(*it) == NOMAD::DS_VAR ) {
      ++it;
      if ( !NOMAD::atoi ( *it , var_index ) ||
	   var_index < 0                    ||
	   var_index >= _dimension             ) {
	return false;
      }
    }
  }
  return true;
}

/*--------------------------------------------*/
/*  check a directory name (static, private)  */
/*--------------------------------------------*/
bool NOMAD::Parameters::check_directory ( std::string & s )
{
  // step 1: remove spaces at the begining:
  size_t i = 0 , ns = s.size();
  while ( i < ns && s[i] == ' ' )
    ++i;
  std::string ss;
  while ( i < ns )
    ss.push_back(s[i++]);
  if ( ss.empty() )
    return false;
  s = ss;

  // step 2: replace '/' or '\' with DIR_SEP:
  i  = 0;
  ns = s.size();
  while ( i < ns ) {
    if ( s[i] == '/' || s[i] == '\\' )
      s[i] = NOMAD::DIR_SEP;
    ++i;
  }

  // step 3: add DIR_SEP at the end:
  if ( i >= 1 ) {
    if (s[i-1] != NOMAD::DIR_SEP) {
      s += NOMAD::DIR_SEP;
    }
  }
  else {
    s = ".";
    s.push_back ( NOMAD::DIR_SEP );
  }

  return true;
}

/*---------------------------------------------------------------*/
/*  interpretation of the Parameter_Entry for PERIODIC_VARIABLE  */
/*  (private)                                                    */
/*---------------------------------------------------------------*/
void NOMAD::Parameters::interpret_periodic_var
( const NOMAD::Parameter_Entries & entries )
{
  int                                      i , j , k;
  std::list<std::string>::const_iterator   it , end;
  NOMAD::Parameter_Entry                 * pe = entries.find ( "PERIODIC_VARIABLE" );

  while ( pe ) {      

    // just one variable index (can be '*' or a range of indexes 'i-j'):
    if ( pe->get_nb_values() == 1 ) {

      it = pe->get_values().begin();
      if ( !NOMAD::string_to_index_range ( *it , i , j , &_dimension ) )
 	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: PERIODIC_VARIABLE" );

      for ( k = i ; k <= j ; ++k )
 	set_PERIODIC_VARIABLE (k);
    }

    // list of variable indexes:
    else {
      end = pe->get_values().end();
      for ( it = pe->get_values().begin() ; it != end ; ++it ) {
 	if ( !NOMAD::atoi ( *it , i ) )
 	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: PERIODIC_VARIABLE" );
	set_PERIODIC_VARIABLE (i);
      }
    }

    pe->set_has_been_interpreted();
    pe = pe->get_next();
  }
}

/*------------------------------------------------------------*/
/*  interpretation of the Parameter_Entry for VARIABLE_GROUP  */
/*  (private)                                                 */
/*------------------------------------------------------------*/
void NOMAD::Parameters::interpret_var_groups ( const NOMAD::Parameter_Entries & entries )
{
  int                                      i , j , k;
  std::set<int>                            var_indexes;
  std::list<std::string>::const_iterator   it , end;
  NOMAD::Parameter_Entry                 * pe = entries.find ( "VARIABLE_GROUP" );

  while ( pe ) {      

    // just one variable index (can be '*' or a range of indexes 'i-j'):
    if ( pe->get_nb_values() == 1 ) {

      it = pe->get_values().begin();
      if ( !NOMAD::string_to_index_range ( *it , i , j , &_dimension ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: VARIABLE_GROUP" );

      for ( k = j ; k >= i ; --k )
	var_indexes.insert(k);
    }

    // list of variable indexes:
    else {
      end = pe->get_values().end();
      for ( it = pe->get_values().begin() ; it != end ; ++it ) {
	if ( !NOMAD::atoi ( *it , i ) )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"invalid parameter: VARIABLE_GROUP" );
	var_indexes.insert(i);
      }
    }

    set_VARIABLE_GROUP ( var_indexes         ,
			 _direction_types    ,
			 _sec_poll_dir_types ,
			 -1                    );

    var_indexes.clear();

    pe->set_has_been_interpreted();
    pe = pe->get_next();
  }
}

/*------------------------------------------------------*/
/*   interpretation of the Parameter_Entry for bounds,  */
/*   fixed variables, and scaling (BFVS)                */
/*   (param_name = "LOWER_BOUND" or "UPPER_BOUND"       */
/*                 "FIXED_VARIABLE" or "SCALING"  )     */
/*   (private)                                          */
/*------------------------------------------------------*/
void NOMAD::Parameters::interpret_BFVS ( const NOMAD::Parameter_Entries & entries    ,
					 const std::string              & param_name   )
{
  // param_name == LOWER_BOUND or UPPER_BOUND or FIXED_VARIABLE:
  if ( param_name != "LOWER_BOUND"    &&
       param_name != "UPPER_BOUND"    &&
       param_name != "FIXED_VARIABLE" &&
       param_name != "SCALING"           )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
		       "wrong use of Parameters::interpret_BFVS()" );
  
  NOMAD::Parameter_Entry                 * pe = entries.find ( param_name );
  std::list<std::string>::const_iterator   it;
  int                                      i, j, k;
  NOMAD::Double                            v;
  std::string                              err;
  std::string                              file_name;
  std::ifstream                            fin;

  while ( pe ) {

    // file name or just one index:
    if ( pe->get_nb_values() == 1 ) {

      // special case for FIXED_VARIABLE without value
      // (the value will be taken from x0, if unique):
      if ( isdigit ( (*pe->get_values().begin())[0] ) ||
	   *pe->get_values().begin() == "*"               ) {

	if ( param_name[0] != 'F' ) {
	  err =  "invalid parameter: " + param_name
  	       + " - only one argument, which is not a file name";
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
	}

	if ( _x0s.size() != 1 || !_x0_cache_file.empty() )
  	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
 		"FIXED_VARIABLE with only a variable index and no unique x0" );
	
  	if ( !NOMAD::string_to_index_range ( *pe->get_values().begin() ,
					     i                         ,
					     j                         ,
					     &_dimension                 ) )
  	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: FIXED_VARIABLE" );

  	for ( k = i ; k <= j ; ++k )
  	  set_FIXED_VARIABLE ( k , (*_x0s[0])[k] );
      }


      // file name:
      else {
	
	file_name = _problem_dir + *pe->get_values().begin();

	fin.open ( file_name.c_str() );
	
	if ( fin.fail() ) {
	  err = "invalid parameter: " + param_name +
	        " - could not open file \'" + file_name + "\'";
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
	}

	try {
	  switch ( param_name[0] ) {
	  case 'L':
	    _lb.reset ( _dimension );
	    fin >> _lb;
	    break;
	  case 'U':
	    _ub.reset ( _dimension );
	    fin >> _ub;
	    break;
	  case 'F':
	    _fixed_variables.reset ( _dimension );
	    fin >> _fixed_variables;
	    break;
	  case 'S':
	    _scaling.reset ( _dimension );
	    fin >> _scaling;
	  }
	}
	catch ( NOMAD::Point::Bad_Input & ) {
	  err = "invalid parameter: " + param_name +
	        " - could not read file \'" + file_name  + "\'";
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
	}

	fin.close();
      }
    }

    // vector form: all values on one row:
    else if ( pe->get_nb_values() == _dimension + 2 ) {

      if ( !pe->is_unique() ) {
	err = "invalid parameter: " + param_name +
	      " - has been given in vector form with [] or () and is not unique";
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }

      it = pe->get_values().begin();
      
      if ( *it != "[" && *it != "(" ) {
	err = "invalid parameter: " + param_name +
	      " - error in vector form with () or []";
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }
      
      ++it;
      for ( k = 0 ; k < _dimension ; ++k ) {
	if ( !v.atof(*it) ) {
	  err = "invalid parameter: " + param_name;
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
	}

	++it;
	switch ( param_name[0] ) {
	case 'L': set_LOWER_BOUND    ( k , v );
	  break;
	case 'U': set_UPPER_BOUND    ( k , v );
	  break;
	case 'F': set_FIXED_VARIABLE ( k , v );
	  break;
	case 'S': set_SCALING        ( k , v );
	}
      }

      if ( *it != "]" && *it != ")" ) {
	err = "invalid parameter: " + param_name +
	      " - error in vector form with () or []";
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }
      
    }

    // indexed values:
    else {

      if ( pe->get_nb_values() != 2 ) {
	err = "invalid parameter: " + param_name;
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }

      it = pe->get_values().begin();
      if ( !NOMAD::string_to_index_range ( *it , i , j , &_dimension ) ) {
	err = "invalid parameter: " + param_name;
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }
      ++it;
      if ( !v.atof(*it) ) {
	err = "invalid parameter: " + param_name;
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }

      for ( k = j ; k >= i ; --k )
	switch (param_name[0]) {
	case 'L': set_LOWER_BOUND    ( k, v );
	  break;
	case 'U': set_UPPER_BOUND    ( k, v );
	  break;
	case 'F': set_FIXED_VARIABLE ( k, v );
	  break;
	case 'S': set_SCALING        ( k, v );
	}
    }
    pe->set_has_been_interpreted();
    pe = pe->get_next();
  }
}

/*----------------------------------------------------------------*/
/*  interpretation of the Parameter_Entry for F_TARGET (private)  */
/*----------------------------------------------------------------*/
void NOMAD::Parameters::interpret_f_target ( const NOMAD::Parameter_Entries & entries )
{
  NOMAD::Double                            d;
  std::list<std::string>::const_iterator   it;
  NOMAD::Parameter_Entry                 * pe = entries.find ( "F_TARGET" );

  if ( pe ) {

    if ( !pe->is_unique() )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: F_TARGET not unique" );

    it = pe->get_values().begin();

    int nb_values = pe->get_nb_values();

    // just one value: single-objective optimization:
    if ( nb_values == 1 ) {

      if ( !d.atof ( *it ) )
 	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: F_TARGET" );
      set_F_TARGET (d);
    }

    // vector form: multi-objective optimization:
    else {
     
      nb_values -= 2;

      NOMAD::Point f_target ( nb_values );

      if ( *it != "[" && *it != "(" )
   	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: F_TARGET - error in vector form with () or []" );
      
      ++it;
      
      for ( int k = 0 ; k < nb_values ; ++k ) {

   	if ( !d.atof ( *it ) )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: F_TARGET" );
	++it;

	f_target[k] = d;
      }

      if ( *it != "]" && *it != ")" )
   	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: F_TARGET - error in vector form with () or []" );

      set_F_TARGET ( f_target );
    }
    pe->set_has_been_interpreted();
  }
}

/*-------------------------------------------------------------*/
/*  interpretation of the Parameter_Entry for mesh/poll sizes  */
/*  (private)                                                  */
/*-------------------------------------------------------------*/
void NOMAD::Parameters::interpret_mesh_sizes
( const NOMAD::Parameter_Entries & entries    ,
  const std::string              & param_name   )
{
  // param_name == "INITIAL_MESH_SIZE" or "MIN_MESH_SIZE" or "MIN_POLL_SIZE":
  if ( param_name != "INITIAL_MESH_SIZE" &&
       param_name != "MIN_MESH_SIZE"     &&
       param_name != "MIN_POLL_SIZE"        )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
		       "wrong use of Parameters::interpret_mesh_sizes()" );

  int                                      i , j , k;
  NOMAD::Double                            v;
  bool                                     relative;
  std::string                              err;
  std::list<std::string>::const_iterator   it;
  NOMAD::Parameter_Entry                 * pe = entries.find ( param_name );

  while ( pe ) {

    // just one value:
    if ( pe->get_nb_values() == 1 ) {

      if ( !pe->is_unique() ) {
	err = "invalid parameter: " + param_name
	    + " - has been given with just one value and is not unique";
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }

      if ( !v.relative_atof ( *pe->get_values().begin() , relative ) ) {
	err = "invalid parameter: " + param_name;
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }

      if ( param_name[0] == 'I' )
	set_INITIAL_MESH_SIZE ( v , relative );
      else if ( param_name[4] == 'M' )
	set_MIN_MESH_SIZE     ( v , relative );
      else
	set_MIN_POLL_SIZE     ( v , relative );
    }

    // indexed form:
    else if ( pe->get_nb_values() == 2 ) {
      
      it = pe->get_values().begin();
      if ( !NOMAD::string_to_index_range ( *it , i , j , &_dimension ) ) {
	err = "invalid parameter: " + param_name;
 	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }
      ++it;

      if ( !v.relative_atof( *it , relative ) ) {
	err = "invalid parameter: " + param_name;
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }

      for ( k = i ; k <= j ; ++k ) {
	if ( param_name[0] == 'I' )
	  set_INITIAL_MESH_SIZE ( k , v , relative );
	else if ( param_name[4] == 'M' )
	  set_MIN_MESH_SIZE     ( k , v , relative );
	else
	  set_MIN_POLL_SIZE     ( k , v , relative );
      }
    }

    // vector form: all values on one row:
    else if ( pe->get_nb_values() == _dimension + 2 ) {

      if ( !pe->is_unique() ) {
	err = "invalid parameter: " + param_name
	    + " - has been given in vector form with [] or () and is not unique";
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }

      it = pe->get_values().begin();
      
      if ( *it != "[" && *it != "(" ) {
	err = "invalid parameter: " + param_name +
	      " - error in vector form with () or []";
  	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }
      
      ++it;
      for ( k = 0 ; k < _dimension ; ++k ) {
  	if ( !v.relative_atof ( *it , relative ) ) {
	  err = "invalid parameter: " + param_name;
  	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
	}
 	++it;
	if ( param_name[0] == 'I' )
	  set_INITIAL_MESH_SIZE ( k , v , relative );
	else if ( param_name[4] == 'M' )
	  set_MIN_MESH_SIZE     ( k , v , relative );
	else
	  set_MIN_POLL_SIZE     ( k , v , relative );
      }

      if ( *it != "]" && *it != ")" ) {
	err = "invalid parameter: " + param_name +
	      " - error in vector form with () or []";
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
      }
    }

    else {
      err = "invalid parameter: " + param_name;
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
    }

    pe->set_has_been_interpreted();
    pe = pe->get_next();
  }
}

/*---------------------------------------------------------------------------------*/
/*          interpretation of the Parameter_Entry for BB_INPUT_TYPE (private)      */
/*---------------------------------------------------------------------------------*/
/*    BB_INPUT_TYPE [ t1 t2 ... tn ]   # blackbox input types (one type/variable)  */
/* or BB_INPUT_TYPE i   t              # ti in { R , C , B , I }                   */
/*                                     #    or { Real , Cat , Bin , Int }          */
/* or BB_INPUT_TYPE i-j t                                                          */
/*---------------------------------------------------------------------------------*/
void NOMAD::Parameters::interpret_bb_input_type
( const NOMAD::Parameter_Entries & entries )
{
  int                                    i , j , k;
  NOMAD::bb_input_type                   bbit;
  std::list<std::string>::const_iterator it;
  NOMAD::Parameter_Entry               * pe = entries.find ( "BB_INPUT_TYPE" );

  while ( pe ) {

    // indexed form:
    if ( pe->get_nb_values() == 2 ) {
      
      it = pe->get_values().begin();
      if ( !NOMAD::string_to_index_range ( *it , i , j , &_dimension ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: BB_INPUT_TYPE" );
      ++it;
      if ( !NOMAD::string_to_bb_input_type ( *it , bbit ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: BB_INPUT_TYPE" );

      for ( k = i ; k <= j ; ++k )
	set_BB_INPUT_TYPE ( k , bbit );
    }

    // vector form: all values on one row:
    else if ( pe->get_nb_values() == _dimension + 2 ) {

      if ( !pe->is_unique() )
 	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	std::string ( "invalid parameter: BB_INPUT_TYPE " )
	+ " - has been given in vector form with [] or () and is not unique" );

      it = pe->get_values().begin();
      
      if ( *it != "[" && *it != "(" )
 	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	"invalid parameter: BB_INPUT_TYPE - error in vector form with () or []" );
      
      ++it;
      for ( k = 0 ; k < _dimension ; ++k ) {
 	if ( !NOMAD::string_to_bb_input_type ( *it , bbit ) )
 	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: BB_INPUT_TYPE" );
	++it;
	set_BB_INPUT_TYPE ( k , bbit );
      }

      if ( *it != "]" && *it != ")" )
 	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	"invalid parameter: BB_INPUT_TYPE - error in vector form with () ot []" );
    }

    else
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				"invalid parameter: BB_INPUT_TYPE" );

    pe->set_has_been_interpreted();
    pe = pe->get_next();
  }
}

/*------------------------------------------------*/
/*  interpretation of the Parameter_Entry for x0  */
/*  (private)                                     */
/*------------------------------------------------*/
void NOMAD::Parameters::interpret_x0 ( const NOMAD::Parameter_Entries & entries )
{
  NOMAD::Parameter_Entry               * pe = entries.find ( "X0" );
  std::list<std::string>::const_iterator it;
  int                                    i , j , k , l;
  NOMAD::Double                          v;
  NOMAD::Point                           tmp_x0;
  std::vector<int>                       indexes;

  while ( pe ) {

    tmp_x0.reset(_dimension);

    // File name:
    if ( pe->get_nb_values() == 1 )
      set_X0 ( *pe->get_values().begin() );

    // Vector form: all values on one row:
    else if ( pe->get_nb_values() == _dimension + 2 ) {
      
      it = pe->get_values().begin();
      
      if ( *it != "[" && *it != "(" ) {

	// particular case with n=1 and 3 entry values:
	// example: X0 1 0 4.0 (first coordinate of the 2nd x0 point put to 4.0)
	if ( _dimension == 1 ) {
	  
	  it = pe->get_values().begin();
	  
	  if ( !NOMAD::atoi ( *it , l ) )
	    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				      "invalid parameter: X0" );

	  i = static_cast<int> ( indexes.size() );
	  if ( l > i )
	    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				      "invalid parameter: X0" );
	  else if ( l == i ) {
	    l = static_cast<int> ( _x0s.size() );
	    indexes.push_back ( l );
	    set_X0 ( tmp_x0 );
	  }
	  else
	    l = indexes[l];
      
	  ++it;
	  if ( !NOMAD::string_to_index_range ( *it , i , j , &_dimension ) )
	    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				      "invalid parameter: X0" );
	  
	  if ( i != 0 && j != 0 )
	    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				      "invalid parameter: X0" );

	  ++it;
  	  if ( !v.atof(*it) )
	    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				      "invalid parameter: X0" );
	  
	  (*_x0s[l])[0] = v;
	}
	
	else
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"invalid parameter: X0 - error in vector form with () or []" );
      }
      
      else {

	++it;
	for ( k = 0 ; k < _dimension ; ++k ) {
   	  if ( !v.atof(*it) )
	    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				      "invalid parameter: X0" );
	  ++it;
	  tmp_x0[k] = v;
	}
	
	if ( *it != "]" && *it != ")" )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"invalid parameter: X0 - error in vector form with () or []" );
	
	set_X0 ( tmp_x0 );
      }
    }

    // Indexed values without x0 index (must be unique)
    // (example: X0 0-5 1.0):
    else if ( pe->get_nb_values() == 2 ) {

      it = pe->get_values().begin();
      if ( !NOMAD::string_to_index_range ( *it , i , j , &_dimension ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , "invalid parameter: X0" );
      ++it;
      if ( !v.atof(*it) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , "invalid parameter: X0" );
      
      if ( indexes.empty() ) {
	l = static_cast<int> ( _x0s.size() );
	indexes.push_back ( l );
	set_X0 ( tmp_x0 );
      }
      else
	l = indexes[0];
      
      for ( k = j ; k >= i ; --k )
	(*_x0s[l])[k] = v;
    }

    // Indexed values with x0 index
    //  example: X0 0 0-5 1.0 --> first x0 point
    //           X0 1 0-5 2.0 --> 2nd x0 point
    else {
      
      if ( pe->get_nb_values() != 3 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , "invalid parameter: X0" );

      it = pe->get_values().begin();
	
      if ( !NOMAD::atoi ( *it , l ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , "invalid parameter: X0" );

      i = static_cast<int> ( indexes.size() );
      if ( l > i )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , "invalid parameter: X0" );
      else if ( l == i ) {
	l = static_cast<int> ( _x0s.size() );
	indexes.push_back ( l );
	set_X0 ( tmp_x0 );
      }
      else
	l = indexes[l];
      
      ++it;
      if ( !NOMAD::string_to_index_range ( *it , i , j , &_dimension ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , "invalid parameter: X0" );
      
      ++it;
      if ( !v.atof(*it) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , "invalid parameter: X0" );
      
      for ( k = j ; k >= i ; --k )
	(*_x0s[l])[k] = v;
      
    }
    pe->set_has_been_interpreted();
    pe = pe->get_next();
  }
}

/*----------------------------------------*/
/*          read a parameters file        */
/*----------------------------------------*/
void NOMAD::Parameters::read ( const std::string & param_file )
{
  // parameters will have to be checked:
  _to_be_checked = true;

  // PROBLEM_DIR:
  // ------------
  _problem_dir.clear();
  size_t k = param_file.find_last_of ( NOMAD::DIR_SEP );
  if ( k >= 0 && k < param_file.size() )
    _problem_dir = param_file.substr (0,k) + NOMAD::DIR_SEP;
  else
    _problem_dir = std::string(".") + NOMAD::DIR_SEP;

  // open the parameters file:
  std::string   err = "could not open parameters file \'" + param_file + "\'";
  std::ifstream fin;
  if ( NOMAD::check_read_file ( param_file ) ) {
    fin.open ( param_file.c_str() );
    if ( !fin.fail() )
      err.clear();
  }
  if ( !err.empty() ) {
    fin.close();
    throw Wrong_Parameters_File ( "Parameters.cpp" , __LINE__ , err );
  }

  // the set of entries:
  NOMAD::Parameter_Entries entries;

  // the file is read: fill the set 'entries' of Parameter_Entry:
  NOMAD::Parameter_Entry * pe;
  std::string              s;

  while ( fin.good() && !fin.eof() ) {

    s.clear();
    
    getline ( fin , s );

    if ( !fin.fail() && !s.empty() ) {
      pe = new NOMAD::Parameter_Entry ( s );
      if ( pe->is_ok() )
	entries.insert ( pe ); // pe will be deleted by ~Parameter_Entries()
      else {
	if ( ( pe->get_name() != "" && pe->get_nb_values() == 0 ) ||
	     pe->get_name() == "STATS_FILE" ) {
	  err = "invalid parameter: " + pe->get_name();
	  delete pe;
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
	}
	delete pe;
      }
    }
  }

  // the file is closed:
  fin.close();

  // entries display:
#ifdef DEBUG
  if ( NOMAD::Slave::is_master() )
    _out << std::endl
	 << NOMAD::open_block ( "parsing of \'" + param_file + "\'" )
	 << entries
	 << NOMAD::close_block();
#endif

  // interpret and set the entries using SET methods:
  std::list<std::string>::const_iterator it , end;
  int                                    i , j , m;
  NOMAD::Double                          d;

  /*----------------------------------------------*/

  // EPSILON:
  // --------
  {
    pe = entries.find ( "EPSILON" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: EPSILON not unique" );
      if ( pe->get_nb_values() != 1 || !d.atof ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: EPSILON" );
      set_EPSILON (d);
      pe->set_has_been_interpreted();     
    }
  }

  // UNDEF_STR:
  // ----------
  pe = entries.find ( "UNDEF_STR" );
  if ( pe ) {
    if ( !pe->is_unique() )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				"invalid parameter: UNDEF_STR not unique" );
    set_UNDEF_STR ( *(pe->get_values().begin()) );
    pe->set_has_been_interpreted();
  }
  
  // INF_STR:
  // --------
  pe = entries.find ( "INF_STR" );
  if ( pe ) {
    if ( !pe->is_unique() )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				"invalid parameter: INF_STR not unique" );
    set_INF_STR ( *(pe->get_values().begin()) );
    pe->set_has_been_interpreted();
  }

  // MESH_UPDATE_BASIS:
  // ------------------
  {
    pe = entries.find ( "MESH_UPDATE_BASIS" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: MESH_UPDATE_BASIS not unique" );
      if ( pe->get_nb_values() != 1 || !d.atof ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MESH_UPDATE_BASIS" );
      set_MESH_UPDATE_BASIS (d);
      pe->set_has_been_interpreted();     
    }
  }

  // INITIAL_MESH_INDEX:
  // -------------------
  {
    pe = entries.find ( "INITIAL_MESH_INDEX" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: INITIAL_MESH_INDEX not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()), i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: INITIAL_MESH_INDEX" );
      pe->set_has_been_interpreted();
      set_INITIAL_MESH_INDEX (i);
    }
  }

  // MESH_REFINING_EXPONENT:
  // -----------------------
  {
    pe = entries.find ( "MESH_REFINING_EXPONENT" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: MESH_REFINING_EXPONENT not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()), i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MESH_REFINING_EXPONENT" );
      pe->set_has_been_interpreted();
      set_MESH_REFINING_EXPONENT (i);
    }
  }

  // MESH_COARSENING_EXPONENT:
  // -------------------------
  {
    pe = entries.find ( "MESH_COARSENING_EXPONENT" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: MESH_COARSENING_EXPONENT not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()), i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MESH_COARSENING_EXPONENT" );
      pe->set_has_been_interpreted();
      set_MESH_COARSENING_EXPONENT (i);
    }
  }

  // MAX_MESH_INDEX:
  // ---------------
  {    
     pe = entries.find ( "MAX_MESH_INDEX" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_MESH_INDEX not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()), i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_MESH_INDEX" );
      pe->set_has_been_interpreted();
      set_MAX_MESH_INDEX (i);
    }
  }

  // HALTON_SEED:
  // ------------
  {
    pe = entries.find ( "HALTON_SEED" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: HALTON_SEED not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()), i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: HALTON_SEED" );
      pe->set_has_been_interpreted();
      set_HALTON_SEED (i);
    }
  }

  // POINT_DISPLAY_LIMIT:
  // --------------------
  {
    pe = entries.find ( "POINT_DISPLAY_LIMIT" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: POINT_DISPLAY_LIMIT not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()), i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: POINT_DISPLAY_LIMIT" );
      set_POINT_DISPLAY_LIMIT (i);
      pe->set_has_been_interpreted();     
    }
  }

  // DIMENSION:
  // ----------
  {
    pe = entries.find ( "DIMENSION" );
    
    if ( !pe ) {
      if ( !pe && _dimension <= 0 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: DIMENSION not defined" );
    }
    else {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: DIMENSION not unique" );

      int dim;
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()), dim) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: DIMENSION" );
      
      pe->set_has_been_interpreted();
      
      set_DIMENSION ( dim );
    }
  }
  
  // SNAP_TO_BOUNDS:
  // ---------------
  {
    pe = entries.find ( "SNAP_TO_BOUNDS" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SNAP_TO_BOUNDS not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SNAP_TO_BOUNDS" );
      set_SNAP_TO_BOUNDS ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // MULTI-MADS:
  // -----------
  {
    // MULTI_OVERALL_BB_EVAL:
    pe = entries.find ( "MULTI_OVERALL_BB_EVAL" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: MULTI_OVERALL_BB_EVAL not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MULTI_OVERALL_BB_EVAL" );
      pe->set_has_been_interpreted();
      set_MULTI_OVERALL_BB_EVAL (i);
    }

    // MULTI_NB_MADS_RUNS:
    pe = entries.find ( "MULTI_NB_MADS_RUNS" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: MULTI_NB_MADS_RUNS not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MULTI_NB_MADS_RUNS" );
      pe->set_has_been_interpreted();
      set_MULTI_NB_MADS_RUNS (i);
    }

    // MULTI_USE_DELTA_CRIT:
    pe = entries.find ( "MULTI_USE_DELTA_CRIT" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: MULTI_USE_DELTA_CRIT not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MULTI_USE_DELTA_CRIT" );
      pe->set_has_been_interpreted();
      set_MULTI_USE_DELTA_CRIT ( i == 1 );
    }    

    // MULTI_F_BOUNDS (f1_min, f2_min, f2_min, f2_max):
    pe = entries.find ( "MULTI_F_BOUNDS" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: MULTI_F_BOUNDS not unique" );
      if ( pe->get_nb_values() != 4 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MULTI_F_BOUNDS" );
      NOMAD::Point mfb ( 4 );
      it = pe->get_values().begin();
      for ( i = 0 ; i < 4 ; ++i ) {
	if ( !d.atof ( *it ) )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: MULTI_F_BOUNDS" );
	mfb[i] = d;
	++it;
      }
      pe->set_has_been_interpreted();
      set_MULTI_F_BOUNDS ( mfb );
    }

    // MULTI_FORMULATION:
    // ------------------
    {
      pe = entries.find ( "MULTI_FORMULATION" );
      if ( pe ) {
	if ( !pe->is_unique() )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"invalid parameter: MULTI_FORMULATION not unique" );
	NOMAD::multi_formulation_type mft;
	if ( pe->get_nb_values() != 1 ||
	     !NOMAD::string_to_multi_formulation_type
	     ( *(pe->get_values().begin()) , mft )    )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"Invalid parameter: MULTI_FORMULATION_TYPE" );
	pe->set_has_been_interpreted();
	set_MULTI_FORMULATION ( mft );
      }
    }
  }

  // SPECULATIVE_SEARCH:
  // -------------------
  {
    pe = entries.find ( "SPECULATIVE_SEARCH" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: SPECULATIVE_SEARCH not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SPECULATIVE_SEARCH" );
      set_SPECULATIVE_SEARCH ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // VNS_SEARCH:
  // -----------
  {
    pe = entries.find ( "VNS_SEARCH" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: VNS_SEARCH not unique" );
      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: VNS_SEARCH" );

      s = *(pe->get_values().begin());
      i = NOMAD::string_to_bool ( s );
      
      // entered as a real:
      if ( i == -1 || s == "1" ) {
	if ( !d.atof ( s ) )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"invalid parameter: VNS_SEARCH" );
	set_VNS_SEARCH ( d );
      }
      // entered as a boolean:
      else
	set_VNS_SEARCH ( i == 1 );

      pe->set_has_been_interpreted();
    }
  }

  // CACHE_SEARCH:
  // -------------
  {
    pe = entries.find ( "CACHE_SEARCH" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: CACHE_SEARCH not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: CACHE_SEARCH" );
      set_CACHE_SEARCH ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // LH_SEARCH:
  // ----------
  {
    pe = entries.find ( "LH_SEARCH" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: LH_SEARCH not unique" );
      if ( pe->get_nb_values() != 2 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: LH_SEARCH" );
      it = pe->get_values().begin();

      if ( !NOMAD::atoi (*it++ , i) || i < 0 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: LH_SEARCH" );

      if ( !NOMAD::atoi (*it , j) || j < 0 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: LH_SEARCH" );

      set_LH_SEARCH ( i , j );
      pe->set_has_been_interpreted();
    }

    // OPPORTUNISTIC_LH:
    // -----------------
    {
      pe = entries.find ( "OPPORTUNISTIC_LH" );
      if ( pe ) {
	if ( !pe->is_unique() )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"invalid parameter: OPPORTUNISTIC_LH not unique" );
	i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
	if ( pe->get_nb_values() != 1 ||  i == -1 )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: OPPORTUNISTIC_LH" );
	set_OPPORTUNISTIC_LH ( i == 1 );
	pe->set_has_been_interpreted();
      }
    }
  }

  // OPPORTUNISTIC_CACHE_SEARCH:
  // ---------------------------
  {
    pe = entries.find ( "OPPORTUNISTIC_CACHE_SEARCH" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: OPPORTUNISTIC_CACHE_SEARCH not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: OPPORTUNISTIC_CACHE_SEARCH" );
      set_OPPORTUNISTIC_CACHE_SEARCH ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }
  
  // opportunistic strategy:
  // -----------------------
  {
    // OPPORTUNISTIC_EVAL:
    pe = entries.find ( "OPPORTUNISTIC_EVAL" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: OPPORTUNISTIC_EVAL not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: OPPORTUNISTIC_EVAL" );
      set_OPPORTUNISTIC_EVAL ( i == 1 );
      pe->set_has_been_interpreted();
    }

    // OPPORTUNISTIC_MIN_NB_SUCCESS:
    pe = entries.find ( "OPPORTUNISTIC_MIN_NB_SUCCESS" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: OPPORTUNISTIC_MIN_NB_SUCCESS not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: OPPORTUNISTIC_MIN_NB_SUCCESS" );
      pe->set_has_been_interpreted();
      set_OPPORTUNISTIC_MIN_NB_SUCCESS (i);
    }

    // OPPORTUNISTIC_MIN_EVAL:
    pe = entries.find ( "OPPORTUNISTIC_MIN_EVAL" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: OPPORTUNISTIC_MIN_EVAL not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: OPPORTUNISTIC_MIN_EVAL" );
      pe->set_has_been_interpreted();
      set_OPPORTUNISTIC_MIN_EVAL (i);
    }

    // OPPORTUNISTIC_MIN_F_IMPRVMT:
    pe = entries.find ( "OPPORTUNISTIC_MIN_F_IMPRVMT" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: OPPORTUNISTIC_MIN_F_IMPRVMT not unique" );
      if ( pe->get_nb_values() != 1 || !d.atof ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: OPPORTUNISTIC_MIN_F_IMPRVMT" );
      pe->set_has_been_interpreted();
      set_OPPORTUNISTIC_MIN_F_IMPRVMT ( d );
    }

    // OPPORTUNISTIC_LUCKY_EVAL:
    pe = entries.find ( "OPPORTUNISTIC_LUCKY_EVAL" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: OPPORTUNISTIC_LUCKY_EVAL not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: OPPORTUNISTIC_LUCKY_EVAL" );
      set_OPPORTUNISTIC_LUCKY_EVAL ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // Directions (DIRECTION_TYPE and SEC_POLL_DIR_TYPE):
  // --------------------------------------------------
  {
    NOMAD::direction_type dt;

    pe = entries.find ( "DIRECTION_TYPE" );
    while ( pe ) {

      if ( !NOMAD::strings_to_direction_type ( pe->get_values() , dt ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: DIRECTION_TYPE" );
      set_DIRECTION_TYPE ( dt );

      pe->set_has_been_interpreted();
      pe = pe->get_next();
    }

    pe = entries.find ( "SEC_POLL_DIR_TYPE" );
    while ( pe ) {
      if ( !NOMAD::strings_to_direction_type ( pe->get_values() , dt ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SEC_POLL_DIR_TYPE" );
      set_SEC_POLL_DIR_TYPE ( dt );

      pe->set_has_been_interpreted();
      pe = pe->get_next();
    }
  }

  // MAX_ITERATIONS:
  // ---------------
  {
    pe = entries.find ( "MAX_ITERATIONS" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_ITERATIONS not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_ITERATIONS" );
      pe->set_has_been_interpreted();
      set_MAX_ITERATIONS (i);
    }
  }

  // MAX_CACHE_MEMORY:
  // -----------------
  {
    pe = entries.find ( "MAX_CACHE_MEMORY" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_CACHE_MEMORY not unique" );
      if ( pe->get_nb_values() != 1 || !d.atof ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_CACHE_MEMORY" );
      pe->set_has_been_interpreted();
      set_MAX_CACHE_MEMORY (static_cast<float>(d.value()));
    }
  }

  // MAX_EVAL:
  // ---------
  {
    pe = entries.find ( "MAX_EVAL" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_EVAL not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_EVAL" );
      pe->set_has_been_interpreted();
      set_MAX_EVAL (i);
    }
  }

  // MAX_BB_EVAL:
  // ------------
  {
    pe = entries.find ( "MAX_BB_EVAL" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_BB_EVAL not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_BB_EVAL" );
      pe->set_has_been_interpreted();
      set_MAX_BB_EVAL (i);
    }
  }

  // MAX_SIM_BB_EVAL:
  // ----------------
  {
    pe = entries.find ( "MAX_SIM_BB_EVAL" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_SIM_BB_EVAL not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_SIM_BB_EVAL" );
      pe->set_has_been_interpreted();
      set_MAX_SIM_BB_EVAL (i);
    }
  }

  // MAX_SGTE_EVAL:
  // --------------
  {
    pe = entries.find ( "MAX_SGTE_EVAL" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_SGTE_EVAL not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()), i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_SGTE_EVAL" );
      pe->set_has_been_interpreted();
      set_MAX_SGTE_EVAL (i);
    }
  }

  // MAX_TIME:
  // ---------
  {
    pe = entries.find ( "MAX_TIME" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_TIME not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()), i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MAX_TIME" );

      pe->set_has_been_interpreted();
      set_MAX_TIME (i);
    }
  }

  // STAT_SUM_TARGET:
  // ----------------
  {
    pe = entries.find ( "STAT_SUM_TARGET" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: STAT_SUM_TARGET not unique" );
      if ( pe->get_nb_values() != 1 || !d.atof ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: STAT_SUM_TARGET" );
      pe->set_has_been_interpreted();
      set_STAT_SUM_TARGET ( d );
    }
  }

  // L_CURVE_TARGET:
  // ---------------
  {
    pe = entries.find ( "L_CURVE_TARGET" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: L_CURVE_TARGET not unique" );
      if ( pe->get_nb_values() != 1 || !d.atof ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: L_CURVE_TARGET" );
      pe->set_has_been_interpreted();
      set_L_CURVE_TARGET ( d );
    }
  }
  
  // EXTENDED_POLL_TRIGGER:
  // ----------------------
  {
    pe = entries.find ( "EXTENDED_POLL_TRIGGER" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: EXTENDED_POLL_TRIGGER not unique" );

      bool rel;

      if ( pe->get_nb_values() != 1 ||
	   !d.relative_atof ( *(pe->get_values().begin()) , rel ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: EXTENDED_POLL_TRIGGER" );

      pe->set_has_been_interpreted();
      set_EXTENDED_POLL_TRIGGER ( d , rel );
    }
  }

  // EXTENDED_POLL_ENABLED:
  // ----------------------
  {
    pe = entries.find ( "EXTENDED_POLL_ENABLED" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: EXTENDED_POLL_ENABLED not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: EXTENDED_POLL_ENABLED" );
      set_EXTENDED_POLL_ENABLED ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // USER_CALLS_ENABLED:
  // -------------------
  {
    pe = entries.find ( "USER_CALLS_ENABLED" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: USER_CALLS_ENABLED not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: USER_CALLS_ENABLED" );
      set_USER_CALLS_ENABLED ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // ASYNCHRONOUS:
  // -------------
  {
    pe = entries.find ( "ASYNCHRONOUS" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: ASYNCHRONOUS not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: ASYNCHRONOUS" );
      set_ASYNCHRONOUS ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // RHO:
  // ----
  {
    pe = entries.find ( "RHO" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: RHO not unique" );
      if ( pe->get_nb_values() != 1 || !d.atof ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: RHO" );
      pe->set_has_been_interpreted();
      set_RHO(d);
    }
  }

  // H_MIN:
  // ------
  {
    pe = entries.find ( "H_MIN" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: H_MIN not unique" );
      if ( pe->get_nb_values() != 1 || !d.atof ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: H_MIN" );
      pe->set_has_been_interpreted();
      set_H_MIN(d);
    }
  }

  // H_MAX_0:
  // --------
  {
    pe = entries.find ( "H_MAX_0" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: H_MAX_0 not unique" );
      if ( pe->get_nb_values() != 1 || !d.atof ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "Invalid parameter: H_MAX_0" );
      pe->set_has_been_interpreted();
      set_H_MAX_0(d);
    }
  }

  // H_NORM:
  // -------
  {
    pe = entries.find ( "H_NORM" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: H_NORM not unique" );
      NOMAD::hnorm_type hn = NOMAD::L2;
      if ( pe->get_nb_values() != 1 ||
	   !NOMAD::string_to_hnorm_type ( *(pe->get_values().begin()) , hn ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "Invalid parameter: H_NORM" );
      pe->set_has_been_interpreted();
      set_H_NORM ( hn );
    }
  }

  // TMP_DIR:
  // --------
  {
    _tmp_dir.clear();
    pe = entries.find ( "TMP_DIR" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: TMP_DIR not unique" );

      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: TMP_DIR" );
      
      set_TMP_DIR ( *(pe->get_values().begin()) );

      pe->set_has_been_interpreted();
    }
  }

  // ADD_SEED_TO_FILE_NAMES:
  // -----------------------
  {
    pe = entries.find ( "ADD_SEED_TO_FILE_NAMES" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: ADD_SEED_TO_FILE_NAMES not unique" );

      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: ADD_SEED_TO_FILE_NAMES" );

      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: ADD_SEED_TO_FILE_NAMES" );
      set_ADD_SEED_TO_FILE_NAMES ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // SOLUTION_FILE:
  // --------------
  {
    _solution_file.clear();
    pe = entries.find ( "SOLUTION_FILE" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SOLUTION_FILE not unique" );

      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SOLUTION_FILE" );
      set_SOLUTION_FILE ( *(pe->get_values().begin()) );
      pe->set_has_been_interpreted();
    }
  }

  // HISTORY_FILE:
  // -------------
  {
    _history_file.clear();
    pe = entries.find ( "HISTORY_FILE" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: HISTORY_FILE not unique" );

      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: HISTORY_FILE" );      
      set_HISTORY_FILE ( *(pe->get_values().begin()) );
      pe->set_has_been_interpreted();
    }
  }

  // STATS_FILE:
  // -----------
  {
    pe = entries.find ( "STATS_FILE" );
    if ( pe ) {

      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: STATS_FILE not unique" );      

      end = pe->get_values().end();
      it  = pe->get_values().begin();
      std::string file_name = *it;
      ++it;

      std::list<std::string> ls;
      while ( it != end ) {
	ls.push_back(*it);
	++it;
      }

      ls.resize(ls.size()-1);

      set_STATS_FILE ( file_name , ls );
      pe->set_has_been_interpreted();
    }
  }

  // CACHE FILE:
  // -----------
  {
    _cache_file.clear();
    pe = entries.find ( "CACHE_FILE" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: CACHE_FILE not unique" );
      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: CACHE_FILE" );
      set_CACHE_FILE ( *(pe->get_values().begin()) );
      pe->set_has_been_interpreted();
    }
  }

  // SGTE_CACHE FILE:
  // ----------------
  {
    _sgte_cache_file.clear();
    pe = entries.find ( "SGTE_CACHE_FILE" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SGTE_CACHE_FILE not unique" );
      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SGTE_CACHE_FILE" );
      set_SGTE_CACHE_FILE ( *(pe->get_values().begin()) );
      pe->set_has_been_interpreted();
    }
  }


  // CACHE_SAVE_PERIOD:
  // ------------------
  {
    pe = entries.find ( "CACHE_SAVE_PERIOD" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: CACHE_SAVE_PERIOD not unique" );
      if ( pe->get_nb_values() != 1 || !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: CACHE_SAVE_PERIOD" );
      set_CACHE_SAVE_PERIOD (i);
      pe->set_has_been_interpreted();
    }
  }

  // SGTE_COST:
  // ----------
  {
    pe = entries.find ( "SGTE_COST" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SGTE_COST not unique" );
      if ( pe->get_nb_values() != 1 ||
	   !NOMAD::atoi (*(pe->get_values().begin()) , i) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SGTE_COST" );
      set_SGTE_COST (i);
      pe->set_has_been_interpreted();
    }
  }

  // X0:
  // ---
  interpret_x0 ( entries );

  // FIXED_VARIABLE:
  // ---------------
  interpret_BFVS ( entries , "FIXED_VARIABLE");

  // LOWER_BOUND:
  // ------------
  interpret_BFVS ( entries , "LOWER_BOUND");

  // UPPER_BOUND:
  // ------------
  interpret_BFVS ( entries , "UPPER_BOUND");

  // SCALING:
  // --------
  interpret_BFVS ( entries , "SCALING" );

  // BB_INPUT_TYPE:
  // --------------
  interpret_bb_input_type ( entries );

  // F_TARGET:
  // ---------
  interpret_f_target ( entries );

  // STOP_IF_FEASIBLE:
  // -----------------
  pe = entries.find ( "STOP_IF_FEASIBLE" );
  if ( pe ) {
    if ( !pe->is_unique() )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: STOP_IF_FEASIBLE not unique" );
    i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
    if ( pe->get_nb_values() != 1 ||  i == -1 )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: STOP_IF_FEASIBLE" );
    pe->set_has_been_interpreted();
    set_STOP_IF_FEASIBLE ( i == 1 );
  }

  // BB_INPUT_INCLUDE_TAG:
  // ---------------------
  {
    pe = entries.find ( "BB_INPUT_INCLUDE_TAG" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: BB_INPUT_INCLUDE_TAG not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: BB_INPUT_INCLUDE_TAG" );
      set_BB_INPUT_INCLUDE_TAG ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // BB_INPUT_INCLUDE_SEED:
  // ----------------------
  {
    pe = entries.find ( "BB_INPUT_INCLUDE_SEED" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: BB_INPUT_INCLUDE_SEED not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: BB_INPUT_INCLUDE_SEED" );
      set_BB_INPUT_INCLUDE_SEED ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }
  // BB_REDIRECTION:
  // ---------------
  {
    pe = entries.find ( "BB_REDIRECTION" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: BB_REDIRECTION not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: BB_REDIRECTION" );
      set_BB_REDIRECTION ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // INITIAL_MESH_SIZE, MIN_MESH_SIZE, and MIN_POLL_SIZE:
  // ----------------------------------------------------
  interpret_mesh_sizes ( entries , "INITIAL_MESH_SIZE" );
  interpret_mesh_sizes ( entries , "MIN_MESH_SIZE"     );
  interpret_mesh_sizes ( entries , "MIN_POLL_SIZE"     );
  
  // BB_OUTPUT_TYPE:
  // ---------------
  {
    pe = entries.find ( "BB_OUTPUT_TYPE" );

    if ( !pe ) {
      if ( _bb_output_type.empty() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: BB_OUTPUT_TYPE not defined" );
    }
    else {

      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: BB_OUTPUT_TYPE not unique" );

      m = pe->get_nb_values();

      if ( m <= 0 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: BB_OUTPUT_TYPE" );
    
      NOMAD::bb_output_type            cur;
      std::list<NOMAD::bb_output_type> bbot;
      i   = 0;
      end = pe->get_values().end();
      for ( it = pe->get_values().begin() ; it != end ; ++it ) {
	if ( !NOMAD::string_to_bb_output_type ( *it , cur ) ) {
	  err = "invalid parameter: BB_OUTPUT_TYPE (" + pe->get_name();
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
	}
	bbot.push_back (cur);
      }

      set_BB_OUTPUT_TYPE ( bbot );
  
      pe->set_has_been_interpreted();
    }
  }

  // BB_EXE:
  // -------
  {
    pe = entries.find ( "BB_EXE" );
    if ( pe ) {

      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: BB_EXE not unique" );

      m = pe->get_nb_values();

      if ( m == 1 )
	set_BB_EXE ( *pe->get_values().begin() );

      else {

	if ( m != static_cast<int>(_bb_output_type.size()) )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: BB_EXE" );
      
	std::list<std::string> bbexe;
	end = pe->get_values().end();
	for ( it = pe->get_values().begin() ; it != end ; ++it )
	  bbexe.push_back (*it);

	set_BB_EXE ( bbexe );
      }

      pe->set_has_been_interpreted();
    }
  }

  // SGTE_EXE:
  // ---------
  {
    pe = entries.find ( "SGTE_EXE" );
    if ( pe ) {

      std::string bb_exe_name , sgte_name;
      
      if ( pe->get_nb_values() == 1 ) {
	if ( !pe->is_unique() )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	        "invalid parameter: SGTE_EXE (with one arguement) not unique" );
	sgte_name = *pe->get_values().begin();
      }

      else if ( pe->get_nb_values() == 2 ) {
	bb_exe_name = *pe->get_values().begin();
	sgte_name   = *(++pe->get_values().begin());
      }

      else
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SGTE_EXE" );

      set_SGTE_EXE ( bb_exe_name , sgte_name );
      pe->set_has_been_interpreted();
    }
  }

  // SGTE_EVAL_SORT:
  // ---------------
  {
    pe = entries.find ( "SGTE_EVAL_SORT" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SGTE_EVAL_SORT not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SGTE_EVAL_SORT" );
      set_SGTE_EVAL_SORT ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // HAS_SGTE:
  // ---------
  {
    pe = entries.find ( "HAS_SGTE" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: HAS_SGTE not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: HAS_SGTE" );
      set_HAS_SGTE ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // OPT_ONLY_SGTE:
  // --------------
  {
    pe = entries.find ( "OPT_ONLY_SGTE" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: OPT_ONLY_SGTE not unique" );
      i = NOMAD::string_to_bool ( *(pe->get_values().begin() ) );
      if ( pe->get_nb_values() != 1 ||  i == -1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: OPT_ONLY_SGTE" );
      set_OPT_ONLY_SGTE ( i == 1 );
      pe->set_has_been_interpreted();
    }
  }

  // DISPLAY_DEGREE:
  // ---------------
  {
    pe = entries.find ( "DISPLAY_DEGREE" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: DISPLAY_DEGREE not unique" );
      if ( pe->get_nb_values() != 1 ||
	   !set_DISPLAY_DEGREE ( *(pe->get_values().begin()) ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: DISPLAY_DEGREE" );
      pe->set_has_been_interpreted();
    }
  }

  // OPEN_BRACE:
  // -----------
  {
    pe = entries.find ( "OPEN_BRACE" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: OPEN_BRACE not unique" );

      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: OPEN_BRACE" );
      
      set_OPEN_BRACE ( *(pe->get_values().begin()) );

      pe->set_has_been_interpreted();
    }
  }

  // CLOSED_BRACE:
  // -------------
  {
    pe = entries.find ( "CLOSED_BRACE" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: CLOSED_BRACE not unique" );

      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: CLOSED_BRACE" );
      
      set_CLOSED_BRACE ( *(pe->get_values().begin()) );

      pe->set_has_been_interpreted();
    }
  }

  // DISPLAY_STATS:
  {
    pe = entries.find ( "DISPLAY_STATS" );
    if ( pe ) {
      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: DISPLAY_STATS not unique" );      
      std::list<std::string> ls;
      end = pe->get_values().end();
      for ( it = pe->get_values().begin() ; it != end ; ++it )
	ls.push_back ( *it );
      ls.resize ( ls.size()-1 );
      set_DISPLAY_STATS ( ls );
      pe->set_has_been_interpreted();
    }
  }

  // SEED:
  // -----
  {
    pe = entries.find ( "SEED" );

    if ( pe ) {

      if ( !pe->is_unique() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SEED not unique" );
      
      if ( pe->get_nb_values() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SEED" );

      s = *(pe->get_values().begin());
      NOMAD::toupper(s);

      if ( s == "NONE" || s == "DIFF" )
	i = -1;
      else if ( !NOMAD::atoi ( s , i ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: SEED" );
      set_SEED ( i );
      pe->set_has_been_interpreted();
    }
  }

  // VARIABLE_GROUP:
  // ---------------
  interpret_var_groups ( entries );
  
  // PERIODIC_VARIABLE:
  // ------------------
  interpret_periodic_var ( entries );

  /*----------------------------------------------*/

  // check the non-interpreted parameters:
  pe = entries.find_non_interpreted();
  if ( pe ) {
    err = "invalid parameter: " + pe->get_name() + " - unknown";
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
  }

  // user must check the parameters with Parameters::check()
}

/*---------------------------------------*/
/*                 display               */
/*---------------------------------------*/
void NOMAD::Parameters::display ( const NOMAD::Display & out ) const
{
  std::list<std::string>::const_iterator it;

  if ( _to_be_checked ) {
    out << "parameters not checked" << std::endl;
    return;
  }

  // problem directory:
  if ( !_problem_dir.empty() ) {
    out << "problem directory    : " << _problem_dir << std::endl;
    if ( _tmp_dir != _problem_dir )
      out << "tmp directory        : " << _tmp_dir << std::endl;
  }

  // dimension:
  out << "dimension            : n=" << _dimension << std::endl;

  // bounds:
  if ( _lb.is_defined() ) {
    out << "lower bounds         : ( ";
    _lb.display ( out , " " , 4 , NOMAD::Point::get_display_limit() );
    out << " )" << std::endl;
  }
  if ( _ub.is_defined() ) {
    out << "upper bounds         : ( ";
    _ub.display ( out , " " , 4 , NOMAD::Point::get_display_limit() );
    out << " )" << std::endl;
  }

  // scaling:
  if ( _scaling.is_defined() ) {
    out << "scaling              : ( ";
    _scaling.display ( out , " " , 4 , NOMAD::Point::get_display_limit() );
    out << " )" << std::endl;
  }

  // fixed variables:
  if ( _fixed_variables.is_defined() ) {
    out << "fixed variables      : ( ";
    _fixed_variables.display ( out , " " , 4 , NOMAD::Point::get_display_limit() );
    out << " )" << std::endl;
  }

  // back-box input types:
  if ( _bb_input_include_tag )
    out << "blackbox input files : include tag" << std::endl;
  if ( _bb_input_include_seed )
    out << "blackbox input files : include seed" << std::endl;

  out << "blackbox input types : ";
  if ( get_signature()->all_continuous() )
    out << "all variables are continuous (R)" << std::endl;
  else
    out << "( " << _bb_input_type << " )" << std::endl;
  
  // extended poll trigger:
  if ( get_signature()->has_categorical() ) {
    if ( _extended_poll_enabled ) {
      out << "extended poll trigger: " << _extended_poll_trigger;
      if ( _relative_ept )
	out << " (relative)";
    }
    else
      out << "extended poll is disabled";
    out << std::endl;
  }
    
  // periodic variables:
  if ( !_periodic_variables.empty() ) {
    out << "periodic variables   : { ";
    for ( size_t k = 0 ; k < _periodic_variables.size() ; ++k )
      if ( _periodic_variables[k] )
	out << k << " ";
    out << "}" << std::endl;
  }
   
  // variable groups:
  if ( _var_groups.size() > 1 ) {
    int i = 0;
    out.open_block ( "variable groups" );
    std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp>::const_iterator
      it2 , end2 = _var_groups.end();
    for ( it2 = _var_groups.begin() ; it2 != end2 ; ++it2 )
      out << NOMAD::open_block ( "group #" + NOMAD::itos ( i++ ) )
	  << **it2 << NOMAD::close_block();
    out.close_block();
  }

  // blackbox outputs:
  {
    bool display_bb_exe   = !_bb_exe.empty();
    bool display_sgte_exe = !_sgte_exe.empty();
    int  m                = static_cast<int>(_bb_output_type.size());
    int  w                = 1+int(log(static_cast<double>(m))/NOMAD::LOG10);
    it = _bb_exe.begin();

    out.open_block ( "blackbox outputs (m=" + NOMAD::itos ( m ) + ")" );
    for ( int i = 0 ; i < m ; ++i ) {
      out << "#" << std::setw(w) << i << " " << std::setw(12) << _bb_output_type[i];
      if ( display_bb_exe ) {
	out << "\t" << *it;
	if ( display_sgte_exe )
	  out << "\t" << get_sgte_exe(*it);
	++it;
      }
      out << std::endl;
    }
    out.close_block();
  }

  // signature (standard or extern):
  out << "signature                       : "
      << ( (_std_signature) ? "standard" : "extern" ) << std::endl;

  // BB_REDIRECTION:
  if ( !_bb_redirection ) {
    out << "blackbox output redirection     : ";
    out.display_yes_or_no ( _bb_redirection );
    out << std::endl;
  }

  // surrogate:
  {
    out << "has surrogate                   : ";
    out.display_yes_or_no ( _has_sgte );
    out << std::endl;
    if ( _has_sgte ) {
      
      // OPT_ONLY_SGTE:
      if ( _opt_only_sgte ) {
	out << "minimize only with surrogate    : ";
	out.display_yes_or_no ( _opt_only_sgte );
	out << std::endl;
      }

      // SGTE_EVAL_SORT:
      out << "sort trial points with surrogate: ";
      out.display_yes_or_no ( _sgte_eval_sort );
      out << std::endl;

      // SGTE_COST:
      out << "surrogate cost                  : ";
      if ( _sgte_cost > 0 )
	out << _sgte_cost
	    << " surrogate evaluations count as one bb evaluation" << std::endl;
      else
	out << "none" << std::endl;
    }
  }
  
  // MULTI-MADS:
  if ( get_nb_obj() > 1 ) {
    out << "multi-MADS                      : [overall bb eval=";
    if ( _multi_overall_bb_eval >= 0 )
      out << _multi_overall_bb_eval;
    else
      out << "-";
    out << "] [nb MADS runs=";
    if ( _multi_nb_mads_runs >= 0 )
      out << _multi_nb_mads_runs;
    else
      out << "-";
    out << "] [use delta crit=";
    out.display_yes_or_no ( _multi_use_delta_crit );
    out << "]" << std::endl
	<< "                                  [formulation="
	<< _multi_formulation << "]";
    if ( _multi_f_bounds.is_defined() ) {
      out << " [f_bounds=";
      _multi_f_bounds.display ( out , "," , -1 , -1 );
      out << "]";
    }
    out << std::endl;
  }

  // barrier:
  if ( _has_constraints ) {

    out << "barrier type                    : ";
    switch ( _barrier_type ) {
      case NOMAD::EB:
      out << "extreme" << std::endl;
      break;
    case NOMAD::PEB_P:
    case NOMAD::PB:
      out << "progressive" << std::endl;
      out << "prog. barrier trigger           : " << _rho << std::endl;
      break;
    default:
      out << "filter" << std::endl;
    }
    out << "barrier h_min                   : " << _h_min    << std::endl
	<< "barrier initial h_max           : " << _h_max_0  << std::endl;
  }
  if ( _has_filter_constraints )
    out << "barrier h_norm                  : " << _h_norm << std::endl;

  // ADD_SEED_TO_FILE_NAMES:
  out << "add seed to output file names   : ";
  out.display_yes_or_no ( _add_seed_to_file_names );
  out << std::endl;

  // SOLUTION_FILE:
  out << "solution file                   : ";
  if ( !_solution_file.empty() )
    out << _solution_file << std::endl;
  else
    out << "none" << std::endl;

  // HISTORY_FILE:
  out << "history file                    : ";
  if ( !_history_file.empty() )
    out << _history_file << std::endl;
  else
    out << "none" << std::endl;

  // STATS_FILE:
  out << "stats file                      : ";
  if ( !_stats_file_name.empty() ) {
    out << "(" << _stats_file_name << ") ";
    std::list<std::string>::const_iterator end = _stats_file.end();
    for ( it = _stats_file.begin() ; it != end ; ++it ) {
      if ( it->empty() )
	out << " ";
      else
	out << *it;
    }
    out << std::endl;
  }
  else
    out << "none" << std::endl;

  // CACHE_FILE:
  out << "cache file                      : ";
  if ( !_cache_file.empty() ) {
    out << _cache_file << std::endl;
    out << "cache save period               : ";
    if ( _cache_save_period <= 0 )
      out << "never";
    else if ( _cache_save_period == 1 )
      out << "every iteration";
    else
      out << "every " << _cache_save_period << " iterations";
    out << std::endl;
  }
  else
    out << "none" << std::endl;

  // surrogate cache file:
  if ( !_sgte_cache_file.empty() )
    out << "surrogate cache file            : "
	<< _sgte_cache_file << std::endl;
 
  // X0:
  if ( _x0s.empty() && _x0_cache_file.empty() )
    out << "x0                              : points in \'"
	<< _cache_file << "\'" << std::endl;
  else {
    bool first = true;
    if ( !_x0_cache_file.empty() ) {
      if ( first ) {
	out << "x0                              : ";
	first = false;
      }
      else
	out << "                                : ";
      out << _x0_cache_file;
      if ( _x0_cache_file != _cache_file )
	out << " (read only)";
      out << std::endl;
    }
    if ( !_x0s.empty() ) {
      size_t x0n = _x0s.size();
      for ( size_t k = 0 ; k < x0n ; ++k ) {
	if ( first ) {
	  out << "x0                              : ";
	  first = false;
	}
	else
	  out << "                                : ";
	out << "( ";
	_x0s[k]->display ( out , " " , 4 , NOMAD::Point::get_display_limit() );
	out << " )" << std::endl;
      }
    }
  }

  // directions:
  {
    std::set<NOMAD::direction_type>::const_iterator it , end = _direction_types.end();
    if ( _direction_types.size() == 1 )
      out << "directions         : "
	  << *_direction_types.begin() << std::endl;
    else {
      out << NOMAD::open_block ( "directions" );
      for ( it = _direction_types.begin() ; it != end ; ++it )
 	out << *it << std::endl;
      out.close_block();
    }
    if ( _barrier_type == NOMAD::PB || _barrier_type == NOMAD::PEB_P ) {
      if ( _sec_poll_dir_types.empty() )
	out << "sec. poll dir. type: no secondary poll" << std::endl;
      else {
	if ( _sec_poll_dir_types.size() == 1 )
	  out << "sec. poll dir. type: "
	      << *_sec_poll_dir_types.begin() << std::endl;
	else {
	  end = _sec_poll_dir_types.end();
	  out << NOMAD::open_block ( "sec. poll dir. types" );
	  for ( it = _sec_poll_dir_types.begin() ; it != end ; ++it )
	    out << *it << std::endl;
	  out.close_block();
	}
      }
    }
    if ( _halton_seed >= 0 )
      out << "user halton seed   : " << _halton_seed << std::endl;
  }

  // mesh:
  {
    out << NOMAD::open_block ( "mesh" )
	<< "update basis       : " << std::setw(3) << _mesh_update_basis << std::endl
	<< "coarsening exponent: " << std::setw(3) << _mesh_coarsening_exponent
	<< std::endl
	<< "refining exponent  : " << std::setw(3) << _mesh_refining_exponent
	<< std::endl
	<< "initial mesh index : " << std::setw(3) << _initial_mesh_index << std::endl;
    if ( _max_mesh_index != NOMAD::UNDEFINED_L )
      out << "max mesh index     : " << std::setw(3) << _max_mesh_index << std::endl;
    out << "initial mesh size  : ( ";
    _initial_mesh_size.display ( out , " " , 4 , NOMAD::Point::get_display_limit() );
    out << " )" << std::endl;
    if ( _min_mesh_size.is_defined() ) {
      out << "min mesh size      : ( ";
      _min_mesh_size.display ( out , " " , 4 , NOMAD::Point::get_display_limit() );
      out << " )" << std::endl;
    }
    if ( _min_poll_size.is_defined() ) {
      out << "min poll size      : ( ";
      _min_poll_size.display ( out , " " , 4 , NOMAD::Point::get_display_limit() );
      out << " )" << std::endl;
    }
    out.close_block();
  }

  // ASYNCHRONOUS:
#ifdef USE_MPI
  out << "asynchronous                     : ";
  out.display_yes_or_no ( _asynchronous );
  out << std::endl;
#endif

  // USER_CALLS_ENABLED:
  if ( !_user_calls_enabled )
    out << "user calls                       : disabled" << std::endl;

  // opportunistic strategy:
  {
    out << "opportunistic evaluations        : ";
    out.display_yes_or_no ( _opportunistic_eval );
    out << std::endl;
    if ( _opportunistic_eval ) {
      if ( _opportunistic_min_nb_success > 0 )
	out << "opportunistic min nb success     : "
	    << _opportunistic_min_nb_success << std::endl;
      if ( _opportunistic_min_eval > 0 )
	out << "opportunistic min nb eval        : "
	    << _opportunistic_min_eval << std::endl;
      if ( _opportunistic_min_f_imprvmt.is_defined() )
	out << "opportunistic min obj improvement: "
	    << _opportunistic_min_f_imprvmt << "%"
	    << std::endl;
      if ( _opportunistic_lucky_eval )
	out << "opportunistic lucky eval         : "
	    << _opportunistic_lucky_eval << std::endl;
    }
  }

  // SNAP_TO_BOUNDS:
  out << "snap to bounds                   : ";
  out.display_yes_or_no ( _snap_to_bounds );
  out << std::endl;

  // SPECULATIVE_SEARCH:
  out << "speculative search               : ";
  out.display_yes_or_no ( _speculative_search );
  out << std::endl;

  // VNS_SEARCH:
  out << "VNS search                       : ";
  out.display_yes_or_no ( _VNS_search );
  if ( _VNS_search )
    out << " (trigger=" << _VNS_trigger << ")";
  out << std::endl;

  // LH_SEARCH:
  out << "Latin-Hypercube (LH) search      : ";
  if ( _LH_search_p0 > 0 || _LH_search_pi > 0 ) {
    out << "#init:"   << _LH_search_p0
	<< ", #iter:" << _LH_search_pi
	<< ", opport:";
    out.display_yes_or_no ( _opportunistic_LH );
  }
  else
    out.display_yes_or_no ( false );
  out << std::endl;

  // CACHE_SEARCH:
  out << "cache search                     : ";
  if ( _cache_search ) {
    out.display_yes_or_no ( true );
    out << ", opport:";
    out.display_yes_or_no ( _opportunistic_cache_search );
  }
  else
    out.display_yes_or_no ( false );
  out << std::endl;

  // random seed / unique tag / run id:
  out << "random seed / run id             : " << _seed << std::endl;

  // EPSILON:
  out << "epsilon                          : "
      << NOMAD::Double::get_epsilon() << std::endl;

  // UNDEF_STR:
  out << "undefined string                 : "
      << NOMAD::Double::get_undef_str() << std::endl;

  // INF_STR:
  out << "infinity string                  : "
      << NOMAD::Double::get_inf_str() << std::endl;

  // DISPLAY_DEGREEs:
  out << NOMAD::open_block ( "display degrees" )
      << "general  : " << _out.get_gen_dd()    << std::endl
      << "search   : " << _out.get_search_dd() << std::endl
      << "poll     : " << _out.get_poll_dd()   << std::endl
      << "iterative: " << _out.get_iter_dd()   << std::endl
      << NOMAD::close_block();

  // DISPLAY_STATS:
  out << "display stats                : ";
  std::list<std::string>::const_iterator end = _display_stats.end();
  for ( it = _display_stats.begin() ; it != end ; ++it ) {
    if ( it->empty() )
      out << " ";
    else
      out << *it;
  }
  out << std::endl;

  // POINT_DISPLAY_LIMIT:
  out << "point display limit          : ";
  if ( NOMAD::Point::get_display_limit() > 0 )
    out << NOMAD::Point::get_display_limit() << std::endl;
  else
    out << "no limit" << std::endl;

  // MAX_EVAL:
  if ( _max_eval > 0 )
    out << "max eval. (bb+cache)         : " << _max_eval << std::endl;

  // MAX_BB_EVAL:
  if ( _max_bb_eval >= 0 ) {
    out << "max number of blackbox eval. : " << _max_bb_eval;
    if ( _max_bb_eval == 0 )
      out << " (no blackbox eval. allowed)";
    out << std::endl;
  }

  // MAX_SIM_BB_EVAL:
  if ( _max_sim_bb_eval >= 0 )
    out << "max simulated blackbox eval. : " << _max_sim_bb_eval << std::endl;

  // MAX_SGTE_EVAL:
  if ( _sgte_max_eval >= 0 ) {
    out << "max surrogate eval.          : " << _sgte_max_eval;
    if ( _sgte_max_eval == 0 )
      out << " (no surrogate eval. allowed)";
    out << std::endl;
  }

  // MAX_ITERATIONS:
  if ( _max_iterations >= 0 ) {
    out << "max iterations               : " << _max_iterations;
    if ( _max_iterations == 0 )
      out << " (no iterations allowed)";
    out << std::endl;
  }

  // MAX_CACHE_MEMORY:
  if ( _max_cache_memory > 0 )
    out << "max cache memory             : " << _max_cache_memory
	<< " MB" << std::endl;

  // MAX_TIME:
  if ( _max_time > 0 )
    out << "max wall-clock time          : " << _max_time << "s" << std::endl;

  // F_TARGET:
  if ( _f_target.is_defined() ) {
    out << "objective target             : ";
    if ( _f_target.size() > 1 ) {
      out << "( ";
      _f_target.display ( out , " " , 4 , -1 );
      out << " )" << std::endl;
    }
    else
      out << _f_target[0] << std::endl;
  }

  // STAT_SUM_TARGET:
  if ( _stat_sum_target.is_defined() )
    out << "stat sum target              : "
	<< _stat_sum_target << std::endl;

  // L_CURVE_TARGET:
  if ( _L_curve_target.is_defined() )
    out << "L-curve target               : "
	<< _L_curve_target << std::endl;

  // STOP_IF_FEASIBLE:
  if ( _stop_if_feasible ) {
    out << "stop if feasible             : ";
    out.display_yes_or_no ( _stop_if_feasible );
    out << std::endl;
  }
}

/*---------------------------------------*/
/*            reset stats file           */
/*---------------------------------------*/
void NOMAD::Parameters::reset_stats_file ( void )
{
  _stats_file.clear();
  _stats_file_name.clear();
}

/*---------------------------------------*/
/*          reset variable groups        */
/*---------------------------------------*/

// 1/2 (private):
void NOMAD::Parameters::reset_variable_groups
( std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> & vg ) const
{
  std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp>::const_iterator end = vg.end() , it;
  for ( it = vg.begin() ; it != end ; ++it )
    delete *it;
  vg.clear();
}

// 2/2 (public):
void NOMAD::Parameters::reset_variable_groups ( void )
{
  _to_be_checked = true;
  reset_variable_groups ( _var_groups      );
  reset_variable_groups ( _user_var_groups );
}

/*---------------------------------------*/
/*          reset fixed variables        */
/*---------------------------------------*/
void NOMAD::Parameters::reset_fixed_variables ( void )
{
  _to_be_checked = true;
  _fixed_variables.clear();
}

/*---------------------------------------*/
/*         reset periodic variables      */
/*---------------------------------------*/
void NOMAD::Parameters::reset_periodic_variables ( void )
{
  _to_be_checked = true;
  _periodic_variables.clear();
}

/*---------------------------------------*/
/*              reset bounds             */
/*---------------------------------------*/
void NOMAD::Parameters::reset_bounds ( void )
{
  _to_be_checked = true;
  _lb.clear();
  _ub.clear();
}

/*---------------------------------------*/
/*             reset scaling             */
/*---------------------------------------*/
void NOMAD::Parameters::reset_scaling ( void )
{
  _to_be_checked = true;
  _scaling.clear();
}

/*----------------------------------------*/
/*            check the parameters        */
/*----------------------------------------*/
void NOMAD::Parameters::check ( bool remove_history_file  ,
				bool remove_solution_file ,
				bool remove_stats_file      )
{  
  if ( !_to_be_checked )
    return;

  int i;
  
  /*--------------------------------------------------*/
  /*  display degree and NOMAD::Point::display_limit  */
  /*--------------------------------------------------*/
  {

#ifdef USE_MPI
    if ( !NOMAD::Slave::is_master() )
      _out.set_degrees ( NOMAD::NO_DISPLAY );
#endif

#ifdef DEBUG
#ifdef USE_MPI
    if ( NOMAD::Slave::is_master() )
#endif
      _out.set_degrees ( NOMAD::FULL_DISPLAY );
#endif

    if ( _out.get_gen_dd() == NOMAD::FULL_DISPLAY )
      set_POINT_DISPLAY_LIMIT ( -1 );
  }

  /*----------------------------*/
  /*          DIMENSION         */
  /*----------------------------*/
  if ( _dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: DIMENSION" );

  /*----------------------------*/
  /*        BB_INPUT_TYPE       */
  /*----------------------------*/
  if ( static_cast<int>(_bb_input_type.size()) != _dimension )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: BB_INPUT_TYPE" );

  /*----------------------------*/
  /*           BOUNDS           */
  /*----------------------------*/
  {
    if ( _lb.size() > _dimension )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				"invalid parameter: LOWER_BOUND" );
    if ( _lb.size() < _dimension )
      _lb.resize ( _dimension );
    
    if ( _ub.size() > _dimension )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				"invalid parameter: UPPER_BOUND" );
    if ( _ub.size() < _dimension )
      _ub.resize ( _dimension );
    
    for ( i = 0 ; i < _dimension ; ++i ) {
      if ( _lb[i].is_defined() && _ub[i].is_defined() ) {
	if ( _lb[i] > _ub[i] )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: LOWER_BOUND or UPPER_BOUND" );
	if ( _lb[i] == _ub[i] )
	  set_FIXED_VARIABLE ( i , _lb[i] );
      }
      
      // integer, binary, and categorical variables:
      if ( _bb_input_type[i] != NOMAD::CONTINUOUS ) {
	
	// binary variables:
	if ( _bb_input_type[i] == NOMAD::BINARY ) {
	  _lb[i] = 0.0;
	  _ub[i] = 1.0;
	}

	// integer and categorical variables:
	else {
	  if ( _lb[i].is_defined() )	  
	    _lb[i] = ceil(_lb[i].value());
	  if ( _ub[i].is_defined() )
	    _ub[i] = floor(_ub[i].value());
	}
      }
    }
  }

  /*----------------------------*/
  /*             Mesh           */
  /*----------------------------*/
  {
    // mesh sizes:
    if ( _initial_mesh_size.size() != _dimension )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				"invalid parameter: INITIAL_MESH_SIZE" );

    // initial mesh size:
    // ------------------
    for ( i = 0 ; i < _dimension ; ++i ) {

      // continuous variables:
      // ---------------------
      if ( _bb_input_type[i] == NOMAD::CONTINUOUS ) {

	// default value for initial mesh size
	// (r0.1 if there are bounds, 1.0 otherwise):
	if ( !_initial_mesh_size[i].is_defined() ) {

	  if ( !_lb[i].is_defined() || !_ub[i].is_defined() )
	    _initial_mesh_size[i] = 1.0;
	  else
	    set_INITIAL_MESH_SIZE ( i , 0.1 , true );
	}

	else if ( _initial_mesh_size[i].value() <  NOMAD::Double::get_epsilon() || 
		  _initial_mesh_size[i].value() <= 0.0                             )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: INITIAL_MESH_SIZE" );
      }

      // binary/categorical variables:
      // -----------------------------
      else if ( _bb_input_type[i] == NOMAD::BINARY      ||
		_bb_input_type[i] == NOMAD::CATEGORICAL    ) {
	_initial_mesh_size[i] = 1.0;
      }
      
      // integer variables:
      // ------------------
      else {

	if ( _initial_mesh_size[i].is_defined() ) {
	  if ( _initial_mesh_size[i] < 1.0 )
	    _initial_mesh_size[i] = 1.0;
	}
	else {
	
	  // default value for initial mesh size
	  // (r0.1 if there are bounds, 1.0 otherwise):
	  if ( !_lb[i].is_defined() || !_ub[i].is_defined() )
	    _initial_mesh_size[i] = 1.0;
	  else
	    set_INITIAL_MESH_SIZE ( i , 0.1 , true );
	}
      }
    }

    // min mesh size \Delta^m_min:
    if ( _min_mesh_size.is_defined() ) {

      if ( _min_mesh_size.size() != _dimension )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MIN_MESH_SIZE" );
      
      for ( i = 0 ; i < _dimension ; ++i )
	if ( _min_mesh_size[i].is_defined() &&
	     (_min_mesh_size[i].value() <  NOMAD::Double::get_epsilon() ||
	      _min_mesh_size[i].value() <= 0.0                             ) )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameters: MIN_MESH_SIZE" );
    }

    // min poll size \Delta^p_min:
    if ( _min_poll_size.is_defined() ) {

      _min_poll_size_defined = true;

      if ( _min_poll_size.size() != _dimension )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: MIN_POLL_SIZE" );
      
      for ( i = 0 ; i < _dimension ; ++i ) {

	// continuous variables:
	if ( _bb_input_type[i] == NOMAD::CONTINUOUS ) {
	  if ( _min_poll_size[i].is_defined() &&
	       (_min_poll_size[i].value() <  NOMAD::Double::get_epsilon() ||
		_min_poll_size[i].value() <= 0.0                             ) )
	    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				      "invalid parameters: MIN_POLL_SIZE" );
	}

	// integer and binary variables:
	else if ( _bb_input_type[i] != NOMAD::CATEGORICAL ) {
	  if ( _min_poll_size[i].is_defined() ) {   
	    if ( _min_poll_size[i] < 1.0 )
	      _min_poll_size[i] = 1.0;
	  }
	  else
	    _min_poll_size[i] = 1.0;
	}
      }
    }
    
    // default min poll size for non-continuous variables:
    else {

      _min_poll_size_defined = false;

      _min_poll_size = NOMAD::Point ( _dimension );
      for ( i = 0 ; i < _dimension ; ++i )
	if ( _bb_input_type[i] != NOMAD::CONTINUOUS  &&
	     _bb_input_type[i] != NOMAD::CATEGORICAL    )
	  _min_poll_size[i] = 1.0;
    }

    // default value for _mesh_update_basis (tau):
    if ( !_mesh_update_basis.is_defined() )
      _mesh_update_basis = 4.0;
    else if ( _mesh_update_basis <= 0.0 )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				"invalid parameters: MESH_UPDATE_BASIS" );

    // compare l0 and lmax:
    if ( _max_mesh_index != NOMAD::UNDEFINED_L && _initial_mesh_index > _max_mesh_index )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameters: MAX_MESH_INDEX or INITIAL_MESH_INDEX" );
  }

  int nb_obj = static_cast<int>(_index_obj.size());

  /*----------------------------*/
  /*         DISPLAY_STATS      */
  /*----------------------------*/
  if ( _display_stats.empty() ) {
    std::list<std::string> ls;
    if ( nb_obj == 1 ) {
      ls.push_back ( NOMAD::Display::get_display_stats_keyword ( NOMAD::DS_BBE ) );
      ls.push_back ( std::string() );
    }
    ls.push_back ( NOMAD::Display::get_display_stats_keyword ( NOMAD::DS_OBJ ) );
    set_DISPLAY_STATS ( ls );
  }

  else if ( !check_display_stats ( _display_stats ) )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: DISPLAY_STATS" );

  /*----------------------------*/
  /*          STATS_FILE        */
  /*----------------------------*/
  if ( !_stats_file_name.empty() ) {
    if ( _stats_file.empty() ) {
      std::list<std::string> ls;
      ls.push_back ( NOMAD::Display::get_display_stats_keyword ( NOMAD::DS_BBE ) );
      ls.push_back ( std::string() );
      ls.push_back ( NOMAD::Display::get_display_stats_keyword ( NOMAD::DS_OBJ ) );
      set_STATS_FILE ( _stats_file_name , ls );
    }
    else if ( !check_display_stats ( _stats_file ) )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				"invalid parameter: STATS_FILE" );
  }

  /*----------------------------*/
  /*           SCALING          */
  /*----------------------------*/
  if ( _scaling.size() > _dimension )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: SCALING" );

  if ( _scaling.size() < _dimension )
    _scaling.resize ( _dimension );

  for ( i = 0; i < _dimension; ++i )
    if ( _scaling[i].is_defined() && _scaling[i] == 0.0 )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: SCALING (zero value)" );

  /*----------------------------*/
  /*       FIXED_VARIABLES      */
  /*----------------------------*/
  if ( _fixed_variables.size() > _dimension )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: FIXED_VARIABLE" );

  if ( _fixed_variables.size() < _dimension )
    _fixed_variables.resize ( _dimension );

  int nb_fixed = 0;
  for ( i = 0; i < _dimension; ++i )
    if ( _fixed_variables[i].is_defined() ) {
      ++nb_fixed;
      if ( (_lb[i].is_defined() && _fixed_variables[i] < _lb[i]) ||
	   (_ub[i].is_defined() && _fixed_variables[i] > _ub[i]) ||
	   ( (_bb_input_type[i] == NOMAD::INTEGER     ||
	      _bb_input_type[i] == NOMAD::CATEGORICAL    )
	     && !_fixed_variables[i].is_integer() )              ||
	   ( _bb_input_type[i] == NOMAD::BINARY && !_fixed_variables[i].is_binary() ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: FIXED_VARIABLE" );
    }

  if ( nb_fixed == _dimension )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: FIXED_VARIABLE - all variables are fixed" );

  _nb_free_variables = _dimension - nb_fixed;

  /*---------------------------*/
  /*      blackbox outputs     */
  /*---------------------------*/
  if ( _bb_output_type.empty() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: BB_OUTPUT_TYPE" );
  if ( _bb_output_type.empty() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: BB_OUTPUT_TYPE - undefined" );
  
  size_t m = _bb_output_type.size();
  
  if ( !_bb_exe.empty() && m != _bb_exe.size() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: BB_EXE: wrong number of blackbox executable names" );

  // surrogate:
  if ( !_sgte_exe.empty() ) {

    _has_sgte = true;

    if ( _bb_exe.empty() )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
      "invalid parameter: SGTE_EXE - no BB_EXE is defined" );

    std::map<std::string,std::string>::const_iterator it;
    std::map<std::string,std::string>::const_iterator end = _sgte_exe.end();
    std::list<std::string>::const_iterator            bb_exe_begin = _bb_exe.begin();
    std::list<std::string>::const_iterator            bb_exe_end = _bb_exe.end();

    // an empty string in _sgte_exe means that there is a unique
    // blackbox with the associated surrogate
    // (SGTE_EXE parameter with only one argument):
    it = _sgte_exe.find("");
    if ( it != end ) {
      if ( _sgte_exe.size() != 1 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
        "invalid parameter: SGTE_EXE - impossible to interpret with one argument" );

      std::string bb_exe_name = *bb_exe_begin;
      std::list<std::string>::const_iterator it2 = ++bb_exe_begin;
      while ( it2 != bb_exe_end ) {
	if ( *it2 != bb_exe_name )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: unique SGTE_EXE without unique blackbox executable" );
	++it2;
      }

      std::string sgte_name = it->second;

      _sgte_exe.clear();
      _sgte_exe[bb_exe_name] = sgte_name;
    }
    else
      for ( it = _sgte_exe.begin() ; it != end ; ++it )
	if ( find ( bb_exe_begin , bb_exe_end , it->first ) == bb_exe_end )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				    "invalid parameter: SGTE_EXE" );
  }
  else if ( !_has_sgte ) {
    _sgte_eval_sort = false;
    _sgte_cost      = -1;
    
    if ( _opt_only_sgte )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				"invalid parameter: OPT_ONLY_SGTE" );
  }

  if ( _opt_only_sgte )
    _sgte_eval_sort = false;

  size_t k;
  
  // _STAT_SUM_ and _STAT_AVG_ checks (each one have to be unique):
  _index_stat_sum = _index_stat_avg = -1;
  for ( k = 0 ; k < m ; ++k ) {
    if ( _bb_output_type[k] == NOMAD::STAT_SUM ) {
      if ( _index_stat_sum >= 0 ) {
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: BB_EXE: more than one STAT_SUM output" );
      }
      _index_stat_sum = static_cast<int>(k);
    }
    else if ( _bb_output_type[k] == NOMAD::STAT_AVG ) {
      if ( _index_stat_avg >= 0 ) {
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: BB_EXE: more than one STAT_AVG output" );
      }
      _index_stat_avg = static_cast<int>(k);
    }
    else if ( _bb_output_type[k] == NOMAD::CNT_EVAL ) {
      if ( _index_cnt_eval >= 0 ) {
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: BB_EXE: more than one CNT_EVAL output" );
      }
      _index_cnt_eval = static_cast<int>(k);
    }
  }

  // F_TARGET:
  if ( _f_target.is_defined() && nb_obj != _f_target.size() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: F_TARGET of bad dimension" );

  /*----------------------------*/
  /*          directions        */
  /*----------------------------*/
  bool use_ortho_mads = false;

  {
    bool use_mads = false;
    std::set<NOMAD::direction_type>::const_iterator it , end = _direction_types.end();

    // default value for primary poll directions:
    if ( _direction_types.empty() ) {
      set_DIRECTION_TYPE ( NOMAD::ORTHO_2N );
      use_mads       = true;
      use_ortho_mads = true;
    }
    else
      for ( it = _direction_types.begin() ; it != end ; ++it ) {
	if ( NOMAD::dir_is_mads ( *it ) )
	  use_mads = true;
	if ( NOMAD::dir_is_orthomads ( *it ) )
	  use_ortho_mads = true;
      }
    
    // default value for secondary poll directions:
    if ( _barrier_type == NOMAD::PB || _barrier_type == NOMAD::PEB_P ) {

      if ( _sec_poll_dir_types.empty() ) {
	if ( use_mads ) {
	  if ( _direction_types.size() == 1 ) {
	    NOMAD::direction_type dt = *(_direction_types.begin());
	    if ( dt == NOMAD::ORTHO_1 || dt == NOMAD::ORTHO_2 )
	      set_SEC_POLL_DIR_TYPE ( NOMAD::ORTHO_1 );
	    else if ( dt == NOMAD::LT_1 || dt == NOMAD::LT_2 )
	      set_SEC_POLL_DIR_TYPE ( NOMAD::LT_1 );
	    else
	      set_SEC_POLL_DIR_TYPE ( (use_ortho_mads) ? NOMAD::ORTHO_2 : NOMAD::LT_2 );
	  }
	  else
	    set_SEC_POLL_DIR_TYPE ( (use_ortho_mads) ? NOMAD::ORTHO_2 : NOMAD::LT_2 );
	}
	else
	  set_SEC_POLL_DIR_TYPE ( NOMAD::GPS_NP1_STATIC );
      }

      else {
	bool old_uom = use_ortho_mads;
	bool old_um  = use_mads;
	end = _sec_poll_dir_types.end();
	for ( it = _sec_poll_dir_types.begin() ; it != end ; ++it ) {
	  if ( *it == NOMAD::NO_DIRECTION ) {
	    _sec_poll_dir_types.clear();
	    use_ortho_mads = old_uom;
	    use_mads       = old_um;
	    break;
	  }
	  if ( NOMAD::dir_is_orthomads (*it) )
	    use_ortho_mads = true;
	  if ( NOMAD::dir_is_mads ( *it ) )
	    use_mads = true;
	}
      }
    }
    else
      _sec_poll_dir_types.clear();

    /*----------------------------*/
    /*     SPECULATIVE_SEARCH     */
    /*----------------------------*/
    if ( !use_mads )
      _speculative_search = false;
  }

  /*----------------------------*/
  /*      periodic variables    */
  /*----------------------------*/
  if ( !_periodic_variables.empty() ) {

    // check the size:
    if ( _dimension != static_cast<int>(_periodic_variables.size()) )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: PERIODIC_VARIABLE - bad size" );

    // check the bounds:
    for ( int k = 0 ; k < _dimension ; ++k )
      if ( _periodic_variables[k] ) {
	if ( !_lb[k].is_defined() )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"invalid parameter: PERIODIC_VARIABLE - lower bound not defined" );
	if ( !_ub[k].is_defined() )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"invalid parameter: PERIODIC_VARIABLE - upper bound not defined" );
      }
  }

  /*----------------------------*/
  /*        variable groups     */
  /*       (and Halton seed)    */
  /*----------------------------*/
  {
    // reset variable groups:
    reset_variable_groups ( _var_groups );

    // Halton seed:
    if ( _halton_seed < 0 ) {
      int * primes = new int [_dimension];
      NOMAD::construct_primes ( _dimension , primes );
      _halton_seed = primes[_dimension-1];
      delete [] primes;
    }

    int def_halton_seed = _halton_seed , halton_seed;
    
    std::vector<bool> in_group ( _dimension );
    for ( i = 0 ; i < _dimension ; ++i )
      in_group[i] = false;

    NOMAD::Variable_Group           * vg;
    std::set<NOMAD::direction_type>   direction_types , sec_poll_dir_types;

    // 1. user groups:
    std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp>::const_iterator
      end = _user_var_groups.end() , it;
    for ( it = _user_var_groups.begin() ; it != end ; ++it ) {

      if ( !(*it)->check ( _fixed_variables , _bb_input_type , &in_group ) )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
				  "invalid parameter: VARIABLE_GROUP" );
      
      halton_seed = (*it)->get_halton_seed();
      if ( use_ortho_mads && halton_seed < 0 )
	halton_seed = def_halton_seed++;

      direction_types    = (*it)->get_direction_types();
      sec_poll_dir_types = (*it)->get_sec_poll_dir_types();

      if ( direction_types.empty() )
	direction_types = _direction_types;

      if ( sec_poll_dir_types.empty() )
	sec_poll_dir_types = _sec_poll_dir_types;

      vg = new NOMAD::Variable_Group ( (*it)->get_var_indexes() ,
				       direction_types          ,
				       sec_poll_dir_types       ,
				       halton_seed              ,
				       _out                       );

      _var_groups.insert ( vg );
    }
    
    // 2. 'automatic' groups for other variables:
    std::set<int> vi_cbi;  // list of cont./bin./int. variables
    std::set<int> vi_cat;  // list of categorical variables

    for ( i = 0 ; i < _dimension ; ++i ) {
      if ( !in_group[i] && !_fixed_variables[i].is_defined() ) {
	if (  _bb_input_type[i] != NOMAD::CATEGORICAL )
	  vi_cbi.insert(i);
	else
	  vi_cat.insert(i);
      }
    }

    // creation of a group for cont./bin./int. variables:
    if ( !vi_cbi.empty() ) {

      halton_seed = -1;
      if ( use_ortho_mads )
	halton_seed = def_halton_seed++; 

      vg = new NOMAD::Variable_Group ( vi_cbi              ,
				       _direction_types    ,
				       _sec_poll_dir_types ,
				       halton_seed         ,
				       _out                  );
      
      _var_groups.insert ( vg );
    }

    // creation of a group for categorical variables:
    if ( !vi_cat.empty() ) {
      vg = new NOMAD::Variable_Group ( vi_cat              ,
				       _direction_types    ,
				       _sec_poll_dir_types ,
				       -1                  ,    // no use of Halton seed
				       _out                  );
      _var_groups.insert ( vg );
    }
  }

  /*----------------------------*/
  /*           TMP_DIR          */
  /*----------------------------*/
  {
    if ( _tmp_dir.empty() )
      _tmp_dir = _problem_dir;

    // check the directory:
    if ( !_tmp_dir.empty() && !NOMAD::check_read_file ( _tmp_dir ) ) {
      std::string err = "invalid parameter: TMP_DIR: cannot access \'" + _tmp_dir + "\'";
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
    }
  }

  /*------------------------------------------------------*/
  /*        SOLUTION_FILE, HISTORY_FILE and STATS_FILE    */
  /*  (depending on the value of ADD_SEED_TO_FILE_NAMES,  */
  /*   the seed is added the file names)                  */
  /*------------------------------------------------------*/
  if ( _add_seed_to_file_names ) {

    std::string s_seed = NOMAD::itos(_seed);
    int         n_seed = static_cast<int>(s_seed.size());

    add_seed_to_file_name ( n_seed , s_seed , _solution_file   );
    add_seed_to_file_name ( n_seed , s_seed , _history_file    );
    add_seed_to_file_name ( n_seed , s_seed , _stats_file_name );
  }

  // remove old history, solution, and stats files:
  std::string old_file;
  if ( remove_history_file && !_history_file.empty() ) {
    old_file = _problem_dir + _history_file;
    remove ( old_file.c_str() );
  }
  if ( remove_stats_file && !_stats_file_name.empty() ) {
    old_file = _problem_dir + _stats_file_name;
    remove ( old_file.c_str() );
  }
  if ( remove_solution_file && !_solution_file.empty() ) {
    old_file = _problem_dir + _solution_file;
    remove ( old_file.c_str() );
  }

  /*----------------------------*/
  /*   opportunistic strategy   */
  /*----------------------------*/
  if ( !_opportunistic_eval ) {
    _sgte_eval_sort               = false;
    _opportunistic_lucky_eval     = false;
    _opportunistic_min_nb_success = -1;
    _opportunistic_min_eval       = -1;
    _opportunistic_min_f_imprvmt.clear();
  }  

  // opportunistic default strategy for LH search:
  //   single-objective: the default is taken the same as OPPORTUNISTIC_EVAL
  //   multi-objective : the default is 'no'
  if ( !_opp_LH_is_defined )
    _opportunistic_LH = ( nb_obj > 1 ) ? false : _opportunistic_eval;

  // opportunistic default strategy for cache search
  // (the same as OPPORTUNISTIC_EVAL):
  if ( !_opp_CS_is_defined )
    _opportunistic_cache_search = false;

  /*----------------------------*/
  /*         MULTI-MADS         */
  /*----------------------------*/
  if ( nb_obj > 1 ) {

    if ( _multi_formulation == NOMAD::UNDEFINED_FORMULATION )
      _multi_formulation = ( _VNS_search ) ? NOMAD::DIST_L2 : NOMAD::PRODUCT;

    if ( _multi_nb_mads_runs < 0 ) {

      if ( _multi_overall_bb_eval < 0 ) {
	_multi_nb_mads_runs = 30;
	if ( !_max_bbe_decided ) {
	  _max_bb_eval = 25 * _nb_free_variables;
	  
	  if ( _LH_search_p0 < 0 )
	    _LH_search_p0 = _max_bb_eval;
	}
      }
      else if ( !_max_bbe_decided ) {
	_max_bb_eval = static_cast<int>
	  ( ceil ( sqrt ( 1.0 * _nb_free_variables * _multi_overall_bb_eval ) ) );
	
	if ( _LH_search_p0 < 0 )
	  _LH_search_p0 = _max_bb_eval;
      }
    }
    else if ( _multi_overall_bb_eval > 0 && !_max_bbe_decided ) {     
      _max_bb_eval = _multi_overall_bb_eval / _multi_nb_mads_runs;
      if ( _multi_nb_mads_runs * _max_bb_eval < _multi_overall_bb_eval )
	++_max_bb_eval;
    }
  }

  /*----------------------------------*/
  /*  signature (standard or extern)  */
  /*----------------------------------*/
  NOMAD::Signature * new_s = new NOMAD::Signature ( _dimension          ,
						    _bb_input_type      ,
						    _initial_mesh_size  ,
						    _min_mesh_size      ,
						    _min_poll_size      ,
						    _lb                 ,
						    _ub                 ,
						    _scaling            ,
						    _fixed_variables    ,
						    _periodic_variables ,
						    _var_groups           );
  // extern signature:
  if ( _extern_signature ) {
    bool fail = ( *new_s != *_extern_signature );
    delete new_s;
    if ( fail )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "Parameters::check(): incompatible extern signature" );
  }

  // standard signature:
  else {
    
    if ( _std_signature ) {
      delete new_s;
      _std_signature->reset ( _dimension          ,
			      _bb_input_type      ,
			      _initial_mesh_size  ,
			      _min_mesh_size      ,
			      _min_poll_size      ,
			      _lb                 ,
			      _ub                 ,
			      _scaling            ,
			      _fixed_variables    ,
			      _periodic_variables ,
			      _var_groups           );
    }
    else {
      _std_signature = new_s;
      _std_signature->set_std();
    }
  }
  
  /*----------------------------*/
  /*              X0            */
  /*----------------------------*/
  {
    if ( _x0s.empty() && _x0_cache_file.empty() ) {
      if ( _LH_search_p0 <= 0 )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "Parameters::check(): no starting point" );
      else {
	_to_be_checked = false;
	NOMAD::Signature * s = get_signature();
	_to_be_checked = true;
	if ( s->has_categorical() )
	  throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
		"Parameters::check(): no starting point with categorical variables" );
      }
    }

    size_t x0n = _x0s.size();
    for ( size_t k = 0 ; k < x0n ; ++k )
      if ( !_x0s[k]->is_complete() )
	throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	      "invalid parameter: x0 with missing coordinates" );

    // avoid _x0_cache_file == _sgte_cache_file :
    if ( !_opt_only_sgte                    &&
	 !_x0_cache_file.empty()            &&
	 !_sgte_cache_file.empty()          &&
	 _x0_cache_file == _sgte_cache_file    )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: x0 and sgte cache file are the same" );
  }

  /*----------------------*/

  _to_be_checked = false;
}

/*-----------------------------------------------------------------*/
/*  add seed to a file name: file_name.ext --> file_name.seed.ext  */
/*  (static, private)                                              */
/*-----------------------------------------------------------------*/
void NOMAD::Parameters::add_seed_to_file_name ( int                 n_seed    ,
						const std::string & s_seed    ,
						std::string       & file_name   )
{ 
  int n_pn = static_cast<int>(file_name.size());

  if ( n_pn == 0 )
    return;

  int         k   = static_cast<int>(file_name.find_last_of("."));
  std::string ext = "";
  std::string fic = file_name;

  if ( k >= 0 && k < n_pn ) {
    fic  = file_name.substr ( 0 , k      );
    ext  = file_name.substr ( k , n_pn-k );
    n_pn = k;
  }
 
  if ( n_pn <= n_seed+1 ||
       fic.substr ( n_pn-n_seed , n_pn-1 ) != s_seed )
    file_name = fic + "." + s_seed + ext;
}

/*----------------------------------------*/
/*               GET methods              */
/*----------------------------------------*/

// get_signature:
NOMAD::Signature * NOMAD::Parameters::get_signature ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_signature(), Parameters::check() must be invoked" );
  if ( !_std_signature && !_extern_signature )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_signature(), no signature is set" );
  return (_std_signature) ? _std_signature : _extern_signature;
}

// get_dimension:
int NOMAD::Parameters::get_dimension ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_dimension(), Parameters::check() must be invoked" );
  return _dimension;
}

// get_nb_free_variables:
int NOMAD::Parameters::get_nb_free_variables ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_nb_free_variables(), Parameters::check() must be invoked" );
  return _nb_free_variables;
}

// get_halton_seed:
int NOMAD::Parameters::get_halton_seed ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_halton_seed(), Parameters::check() must be invoked" );
  return _halton_seed;
}

// get_add_seed_to_file_names:
bool NOMAD::Parameters::get_add_seed_to_file_names ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_add_seed_to_file_names(), Parameters::check() must be invoked" );
  return _add_seed_to_file_names;
}

// get_snap_to_bounds:
bool NOMAD::Parameters::get_snap_to_bounds ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_snap_to_bounds(), Parameters::check() must be invoked" );
  return _snap_to_bounds;
}

// get_speculative_search:
bool NOMAD::Parameters::get_speculative_search ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_speculative_search(), Parameters::check() must be invoked" );
  return _speculative_search;
}

// get_cache_search:
bool NOMAD::Parameters::get_cache_search ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_cache_search(), Parameters::check() must be invoked" );
  return _cache_search;
}

// get_VNS_search:
bool NOMAD::Parameters::get_VNS_search ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
          "Parameters::get_VNS_search(), Parameters::check() must be invoked" );
  return _VNS_search;
}

// get_VNS_trigger:
const NOMAD::Double & NOMAD::Parameters::get_VNS_trigger ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
          "Parameters::get_VNS_trigger(), Parameters::check() must be invoked" );
  return _VNS_trigger;
}

// get_LH_search_p0:
int NOMAD::Parameters::get_LH_search_p0 ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
          "Parameters::get_LH_search_p0(), Parameters::check() must be invoked" );
  return _LH_search_p0;
}

// get_LH_search_p0:
int NOMAD::Parameters::get_LH_search_pi ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
          "Parameters::get_LH_search_pi(), Parameters::check() must be invoked" );
  return _LH_search_pi;
}

// get_direction_types:
const std::set<NOMAD::direction_type> &
NOMAD::Parameters::get_direction_types ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
          "Parameters::get_direction_types(), Parameters::check() must be invoked" );
  return _direction_types;
}

// get_sec_poll_dir_types:
const std::set<NOMAD::direction_type> &
NOMAD::Parameters::get_sec_poll_dir_types ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_sec_poll_dir_types(), Parameters::check() must be invoked" );
  return _sec_poll_dir_types;
}

// check if there are Ortho-MADS directions:
bool NOMAD::Parameters::has_orthomads_directions ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::has_orthomads_directions(), Parameters::check() must be invoked" );
  bool use_ortho_mads = NOMAD::dirs_are_orthomads ( _direction_types );
  if ( !use_ortho_mads )
    use_ortho_mads = NOMAD::dirs_are_orthomads ( _sec_poll_dir_types );
  return use_ortho_mads;
}

// get_mesh_update_basis:
const NOMAD::Double & NOMAD::Parameters::get_mesh_update_basis ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_mesh_update_basis(), Parameters::check() must be invoked" );
  return _mesh_update_basis;
}

// get_mesh_coarsening_exponent:
int NOMAD::Parameters::get_mesh_coarsening_exponent ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
   "Parameters::get_mesh_coarsening_exponent(), Parameters::check() must be invoked" );
  return _mesh_coarsening_exponent;
}

// get_mesh_refining_exponent:
int NOMAD::Parameters::get_mesh_refining_exponent ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_mesh_refining_exponent(), Parameters::check() must be invoked" );
  return _mesh_refining_exponent;
}

// get_initial_mesh_index:
int NOMAD::Parameters::get_initial_mesh_index ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_initial_mesh_index(), Parameters::check() must be invoked" );
  return _initial_mesh_index;
}

// get_max_mesh_index:
int NOMAD::Parameters::get_max_mesh_index ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
          "Parameters::get_max_mesh_index(), Parameters::check() must be invoked" );
  return _max_mesh_index;
}

// get_initial_mesh_size:
const NOMAD::Point & NOMAD::Parameters::get_initial_mesh_size ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_initial_mesh_size(), Parameters::check() must be invoked" );
  return _initial_mesh_size;
}

// get_min_mesh_size:
const NOMAD::Point & NOMAD::Parameters::get_min_mesh_size ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_min_mesh_size(), Parameters::check() must be invoked" );
  return _min_mesh_size;
}

// get_min_poll_size:
const NOMAD::Point & NOMAD::Parameters::get_min_poll_size ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_min_poll_size(), Parameters::check() must be invoked" );
  return _min_poll_size;
}

// get_min_poll_size_defined:
bool NOMAD::Parameters::get_min_poll_size_defined ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_min_poll_size_defined(), Parameters::check() must be invoked" );
  return _min_poll_size_defined;
}

// get_extended_poll_trigger:
const NOMAD::Double & NOMAD::Parameters::get_extended_poll_trigger ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_extended_poll_trigger(), Parameters::check() must be invoked" );
  return _extended_poll_trigger;
}

// get_relative_ept:
bool NOMAD::Parameters::get_relative_ept ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_relative_ept(), Parameters::check() must be invoked" );
  return _relative_ept;
}

// get_extended_poll_enabled:
bool NOMAD::Parameters::get_extended_poll_enabled ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_extended_poll_enabled(), Parameters::check() must be invoked" );
  return _extended_poll_enabled;
}

// get_user_calls_enabled:
bool NOMAD::Parameters::get_user_calls_enabled ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_user_calls_enabled(), Parameters::check() must be invoked" );
  return _user_calls_enabled;
}

// get_asynchronous:
bool NOMAD::Parameters::get_asynchronous ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_asynchronous(), Parameters::check() must be invoked" );
  return _asynchronous;
}

// get_x0s:
const std::vector<NOMAD::Point *> & NOMAD::Parameters::get_x0s ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
		       "Parameters::get_x0s(), Parameters::check() must be invoked" );
  return _x0s;
}

// get_x0_cache_file:
const std::string & NOMAD::Parameters::get_x0_cache_file ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_x0_cache_file(), Parameters::check() must be invoked" );
  return _x0_cache_file;
}

// get_lb:
const NOMAD::Point & NOMAD::Parameters::get_lb ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
		       "Parameters::get_lb(), Parameters::check() must be invoked" );
  return _lb;
}

// get_ub:
const NOMAD::Point & NOMAD::Parameters::get_ub ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
		       "Parameters::get_ub(), Parameters::check() must be invoked" );
  return _ub;
}

// get_scaling:
const NOMAD::Point & NOMAD::Parameters::get_scaling ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_scaling(), Parameters::check() must be invoked" );
  return _scaling;
}

// get_fixed_variables:
const NOMAD::Point & NOMAD::Parameters::get_fixed_variables ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_fixed_variables(), Parameters::check() must be invoked" );
  return _fixed_variables;
}

// variable_is_fixed:
bool NOMAD::Parameters::variable_is_fixed ( int index ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::variable_is_fixed(), Parameters::check() must be invoked" );
  if ( index < 0 || index >= _fixed_variables.size() )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::variable_is_fixed(), bad variable index" );
  return _fixed_variables[index].is_defined();
}

// get_bb_nb_outputs:
int NOMAD::Parameters::get_bb_nb_outputs ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_bb_nb_outputs(), Parameters::check() must be invoked" );
  return static_cast<int>(_bb_output_type.size());
}

// get_bb_exe:
const std::list<std::string> & NOMAD::Parameters::get_bb_exe ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
          "Parameters::get_bb_exe(), Parameters::check() must be invoked" );
  return _bb_exe;
}

// get_sgte_eval_sort:
bool NOMAD::Parameters::get_sgte_eval_sort ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_sgte_eval_sort(), Parameters::check() must be invoked" );
  return _sgte_eval_sort;
}

// has_sgte:
bool NOMAD::Parameters::has_sgte ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::has_sgte(), Parameters::check() must be invoked" );
  return _has_sgte;
}

// get_opt_only_sgte:
bool NOMAD::Parameters::get_opt_only_sgte ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_opt_only_sgte(), Parameters::check() must be invoked" );
  return _opt_only_sgte;
}

// has_sgte_exe:
bool NOMAD::Parameters::has_sgte_exe ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::has_sgte_exe(), Parameters::check() must be invoked" );
  return !_sgte_exe.empty();
}

// get_sgte_cost:
int NOMAD::Parameters::get_sgte_cost ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_sgte_cost(), Parameters::check() must be invoked" );
  return _sgte_cost;
}

// get_sgte_exe (returns an empty string if bb_exe has no surrogate):
std::string NOMAD::Parameters::get_sgte_exe ( const std::string & bb_exe ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_sgte_exe(), Parameters::check() must be invoked" );
  std::map<std::string,std::string>::const_iterator it = _sgte_exe.find(bb_exe);
  std::string s;
  if ( it != _sgte_exe.end() )
    s = it->second;
  return s;
}

// get_index_obj:
const std::list<int> & NOMAD::Parameters::get_index_obj ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_index_obj(), Parameters::check() must be invoked" );
  return _index_obj;
}

// get_nb_obj:
int NOMAD::Parameters::get_nb_obj ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_nb_obj(), Parameters::check() must be invoked" );
  return static_cast<int>(_index_obj.size());
}

// get_bb_input_include_tag:
bool NOMAD::Parameters::get_bb_input_include_tag ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_bb_input_include_tag(), Parameters::check() must be invoked" );
  return _bb_input_include_tag;
}

// get_bb_input_include_seed:
bool NOMAD::Parameters::get_bb_input_include_seed ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_bb_input_include_seed(), Parameters::check() must be invoked" );
  return _bb_input_include_seed;
}

// get_bb_redirection:
bool NOMAD::Parameters::get_bb_redirection ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_bb_redirection(), Parameters::check() must be invoked" );
  return _bb_redirection;
}

// get_bb_input_type:
const std::vector<NOMAD::bb_input_type> &
NOMAD::Parameters::get_bb_input_type ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_bb_input_type(), Parameters::check() must be invoked" );
  return _bb_input_type;
}

// get_bb_output_type:
const std::vector<NOMAD::bb_output_type> &
NOMAD::Parameters::get_bb_output_type ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_bb_output_type(), Parameters::check() must be invoked" );
  return _bb_output_type;
}

// get_seed:
int NOMAD::Parameters::get_seed ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_seed(), Parameters::check() must be invoked" );
  return _seed;
}

// get_display_stats:
const std::list<std::string> & NOMAD::Parameters::get_display_stats ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_display_stats(), Parameters::check() must be invoked" );
  return _display_stats;
}

// get_stats_file_name:
const std::string & NOMAD::Parameters::get_stats_file_name ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_stats_file_name(), Parameters::check() must be invoked" );
  return _stats_file_name;
}

// get_stats_file:
const std::list<std::string> & NOMAD::Parameters::get_stats_file ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_stats_file(), Parameters::check() must be invoked" );
  return _stats_file;
}

// get_point_display_limit:
int NOMAD::Parameters::get_point_display_limit ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_point_display_limit(), Parameters::check() must be invoked" );
  return NOMAD::Point::get_display_limit();
}

// out (ex get_display()):
NOMAD::Display & NOMAD::Parameters::out ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
		       "Parameters::out(), Parameters::check() must be invoked" );
  return _out;
}

// get_display_degree 1/2:
void NOMAD::Parameters::get_display_degree ( std::string & d ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_display_degree(), Parameters::check() must be invoked" );
  _out.get_display_degree ( d );
}

// get_display_degree 2/2:
int NOMAD::Parameters::get_display_degree ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_display_degree(), Parameters::check() must be invoked" );
  return _out.get_gen_dd();
}

// get_max_eval:
int NOMAD::Parameters::get_max_eval ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_max_eval(), Parameters::check() must be invoked" );
  return _max_eval;
}

// get_max_bb_eval:
int NOMAD::Parameters::get_max_bb_eval ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_max_bb_eval(), Parameters::check() must be invoked" );
  return _max_bb_eval;
}

// get_max_sim_bb_eval:
int NOMAD::Parameters::get_max_sim_bb_eval ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_max_sim_bb_eval(), Parameters::check() must be invoked" );
  return _max_sim_bb_eval;
}

// get_max_sgte_eval:
int NOMAD::Parameters::get_max_sgte_eval ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_max_sgte_eval(), Parameters::check() must be invoked" );
  return _sgte_max_eval;
}

// get_max_time:
int NOMAD::Parameters::get_max_time ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_max_time(), Parameters::check() must be invoked" );
  return _max_time;
}

// get_max_iterations:
int NOMAD::Parameters::get_max_iterations ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_max_iterations(), Parameters::check() must be invoked" );
  return _max_iterations;
}

// get_max_cache_memory:
float NOMAD::Parameters::get_max_cache_memory ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_max_cache_memory(), Parameters::check() must be invoked" );
  return _max_cache_memory;
}

// get_cache_save_period:
int NOMAD::Parameters::get_cache_save_period ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_cache_save_period(), Parameters::check() must be invoked" );
  return _cache_save_period;
}

// get_stop_if_feasible:
bool NOMAD::Parameters::get_stop_if_feasible ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_stop_if_feasible(), Parameters::check() must be invoked" );
  return _stop_if_feasible;
}

// get_f_target:
const NOMAD::Point & NOMAD::Parameters::get_f_target ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_f_target(), Parameters::check() must be invoked" );
  return _f_target;
}

// get_stat_sum_target:
const NOMAD::Double & NOMAD::Parameters::get_stat_sum_target ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_stat_sum_target(), Parameters::check() must be invoked" );
  return _stat_sum_target;
}

// get_L_curve_target:
const NOMAD::Double & NOMAD::Parameters::get_L_curve_target ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_L_curve_target(), Parameters::check() must be invoked" );
  return _L_curve_target;
}

// get_problem_dir:
const std::string & NOMAD::Parameters::get_problem_dir ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_problem_dir(), Parameters::check() must be invoked" );
  return _problem_dir;
}

// get_tmp_dir:
const std::string & NOMAD::Parameters::get_tmp_dir ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_tmp_dir(), Parameters::check() must be invoked" );
  return _tmp_dir;
}

// get_solution_file:
const std::string & NOMAD::Parameters::get_solution_file ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_solution_file(), Parameters::check() must be invoked" );
  return _solution_file;
}

// get_history_file:
const std::string & NOMAD::Parameters::get_history_file ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_history_file(), Parameters::check() must be invoked" );
  return _history_file;
}

// get_cache_file:
const std::string & NOMAD::Parameters::get_cache_file ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_cache_file(), Parameters::check() must be invoked" );
  return _cache_file;
}

// get_sgte_cache_file:
const std::string & NOMAD::Parameters::get_sgte_cache_file ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_sgte_cache_file(), Parameters::check() must be invoked" );
  return _sgte_cache_file;
}

// get_rho:
const NOMAD::Double & NOMAD::Parameters::get_rho ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
		       "Parameters::get_rho(), Parameters::check() must be invoked" );
  return _rho;
}

// get_h_min:
const NOMAD::Double & NOMAD::Parameters::get_h_min ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_h_min(), Parameters::check() must be invoked" );
  return _h_min;
}

// get_h_max_0:
const NOMAD::Double & NOMAD::Parameters::get_h_max_0 ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_h_max_0(), Parameters::check() must be invoked" );
  return _h_max_0;
}

// get_h_norm:
NOMAD::hnorm_type NOMAD::Parameters::get_h_norm ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_h_norm(), Parameters::check() must be invoked" );
  return _h_norm;
}

// use_sec_poll_center:
bool NOMAD::Parameters::use_sec_poll_center ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::use_second_poll_center(), Parameters::check() must be invoked" );
  return _barrier_type == NOMAD::PB || _barrier_type == NOMAD::PEB_P;
}

// get_barrier_type:
NOMAD::bb_output_type NOMAD::Parameters::get_barrier_type ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_filter_type(), Parameters::check() must be invoked" );
  return _barrier_type;
}

// has_constraints:
bool NOMAD::Parameters::has_constraints ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::has_constraints(), Parameters::check() must be invoked" );
  return _has_constraints;
}

// has_EB_constraints:
bool NOMAD::Parameters::has_EB_constraints ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::has_EB_constraints(), Parameters::check() must be invoked" );
  return _has_EB_constraints;
}

// get_multi_nb_mads_runs:
int NOMAD::Parameters::get_multi_nb_mads_runs ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_multi_nb_mads_runs(), Parameters::check() must be invoked" );
  return _multi_nb_mads_runs;
}

// get_multi_overall_bb_eval:
int NOMAD::Parameters::get_multi_overall_bb_eval ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_multi_overall_bb_eval(), Parameters::check() must be invoked" );
  return _multi_overall_bb_eval;
}

// get_multi_use_delta_crit:
bool NOMAD::Parameters::get_multi_use_delta_crit ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_multi_use_delta_crit(), Parameters::check() must be invoked" );
  return _multi_use_delta_crit;
}

// get_multi_f_bounds:
const NOMAD::Point & NOMAD::Parameters::get_multi_f_bounds ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_multi_f_bounds(), Parameters::check() must be invoked" );
  return _multi_f_bounds;
}

// get_multi_formulation:
NOMAD::multi_formulation_type NOMAD::Parameters::get_multi_formulation ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_multi_formulation(), Parameters::check() must be invoked" );
  return _multi_formulation;
}

// get_opportunistic_cache_search:
bool NOMAD::Parameters::get_opportunistic_cache_search ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
  "Parameters::get_opportunistic_cache_search(), Parameters::check() must be invoked" );
  return _opportunistic_cache_search;
}

// get_opportunistic_LH:
bool NOMAD::Parameters::get_opportunistic_LH ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::get_opportunistic_LH(), Parameters::check() must be invoked" );
  return _opportunistic_LH;
}

// get_opportunistic_eval:
bool NOMAD::Parameters::get_opportunistic_eval ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_opportunistic_eval(), Parameters::check() must be invoked" );
  return _opportunistic_eval;
}

// get_opportunistic_min_nb_success
int NOMAD::Parameters::get_opportunistic_min_nb_success ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
  "Parameters::get_opportunistic_min_nb_success(), Parameters::check() must be invoked");
  return _opportunistic_min_nb_success;
}

// get_opportunistic_min_eval
int NOMAD::Parameters::get_opportunistic_min_eval ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_opportunistic_min_eval(), Parameters::check() must be invoked" );
  return _opportunistic_min_eval;
}

// get_opportunistic_min_f_imprvmt
const NOMAD::Double & NOMAD::Parameters::get_opportunistic_min_f_imprvmt ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
  "Parameters::get_opportunistic_min_f_imprvmt(), Parameters::check() must be invoked" );
  return _opportunistic_min_f_imprvmt;
}

// get_opportunistic_lucky_eval:
bool NOMAD::Parameters::get_opportunistic_lucky_eval ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_opportunistic_lucky_eval(), Parameters::check() must be invoked" );
  return _opportunistic_lucky_eval;
}

// check_stat_sum:
bool NOMAD::Parameters::check_stat_sum ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::check_stat_sum(), Parameters::check() must be invoked" );
  return ( _index_stat_sum >= 0 );
}

// check_stat_avg:
bool NOMAD::Parameters::check_stat_avg ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::check_stat_avg(), Parameters::check() must be invoked" );
  return ( _index_stat_avg >= 0 );
}

// get_index_stat_sum:
int NOMAD::Parameters::get_index_stat_sum ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_index_stat_sum(), Parameters::check() must be invoked" );
  return _index_stat_sum;
}

// get_index_stat_avg:
int NOMAD::Parameters::get_index_stat_avg ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_index_stat_avg(), Parameters::check() must be invoked" );
  return _index_stat_avg;
}

// get_index_cnt_eval:
int NOMAD::Parameters::get_index_cnt_eval ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_index_cnt_eval(), Parameters::check() must be invoked" );
  return _index_cnt_eval;
}

// get_periodic_variables:
const std::vector<bool> & NOMAD::Parameters::get_periodic_variables ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_periodic_variables(), Parameters::check() must be invoked" );
  return _periodic_variables;
}

// get_variable_groups:
const std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> &
NOMAD::Parameters::get_variable_groups ( void ) const
{
  if ( _to_be_checked )
    throw Bad_Access ( "Parameters.cpp" , __LINE__ ,
    "Parameters::get_variable_groups(), Parameters::check() must be invoked" );
  return _var_groups;
}

/*----------------------------------------*/
/*               SET methods              */
/*----------------------------------------*/

// set_POINT_DISPLAY_LIMIT:
void NOMAD::Parameters::set_POINT_DISPLAY_LIMIT ( int dl )
{
  NOMAD::Point::set_display_limit ( dl );
}

// set_DIMENSION:
bool NOMAD::Parameters::set_DIMENSION ( int dim )
{
  if ( _dimension > 0 ) {
    _dimension = -1;
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: DIMENSION - defined twice" );
    return false;
  }

  _to_be_checked = true;
  _dimension     = dim;
  if ( _dimension <= 0 ) {
    _dimension = -1;
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: DIMENSION" );
    return false;
  }

  // all variables are initially considered continuous:
  _bb_input_type.resize ( _dimension );
  for ( int i = 0 ; i < _dimension ; ++i )
    _bb_input_type[i] = NOMAD::CONTINUOUS;

  // resize of _initial_mesh_size:
  _initial_mesh_size.reset ( _dimension );

  return true;
}

// set_EXTERN_SIGNATURE (set a new extern signature and
//                       delete the standard signature):
void NOMAD::Parameters::set_EXTERN_SIGNATURE ( NOMAD::Signature * s )
{
  if ( _std_signature && s == _std_signature )
    return;

  // standard signature:
  delete _std_signature;
  _std_signature    = NULL;
  _extern_signature = s;

  // dimension:
  _dimension = -1;
  set_DIMENSION ( s->get_n() );

  // input types:
  set_BB_INPUT_TYPE ( s->get_input_types() );

  // bounds:
  set_LOWER_BOUND ( s->get_lb() );
  set_UPPER_BOUND ( s->get_ub() );
	
  // scaling:
  set_SCALING ( s->get_scaling() );

  // fixed variables:
  set_FIXED_VARIABLE ( s->get_fixed_variables() );

  // periodic variables:
  set_PERIODIC_VARIABLE ( s->get_periodic_variables() );

  // variable groups:
  reset_variable_groups();
  set_VARIABLE_GROUP ( s->get_var_groups() ); 

  _to_be_checked = true;
}

// set_HALTON_SEED:
void NOMAD::Parameters::set_HALTON_SEED ( int hs )
{
  _to_be_checked = true;
  _halton_seed   = ( hs < 0 ) ? -1 : hs;
}

// set_SNAP_TO_BOUNDS:
void NOMAD::Parameters::set_SNAP_TO_BOUNDS ( bool stb )
{
  _to_be_checked  = true;
  _snap_to_bounds = stb;
}

// set_SPECULATIVE_SEARCH:
void NOMAD::Parameters::set_SPECULATIVE_SEARCH ( bool ss )
{
  _to_be_checked      = true;
  _speculative_search = ss;
}

// set_CACHE_SEARCH:
void NOMAD::Parameters::set_CACHE_SEARCH ( bool s )
{
  _to_be_checked = true;
  _cache_search  = s;
}

// set_VNS_SEARCH (1/2):
void NOMAD::Parameters::set_VNS_SEARCH ( bool s )
{
  _to_be_checked = true;
  _VNS_search    = s;
  _VNS_trigger   = ( s ) ? 0.75 : NOMAD::Double();
}

// set_VNS_SEARCH (2/2):
void NOMAD::Parameters::set_VNS_SEARCH ( const NOMAD::Double & trigger )
{
  _to_be_checked = true;
  if ( !trigger.is_defined() ) {
    _VNS_search = false;
    return;
  }

  if ( trigger < 0.0 || trigger > 1.0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: VNS_SEARCH: must be in [0;1]" );

  _VNS_search  = ( trigger > 0.0 );
  _VNS_trigger = trigger;
}

// set_LH_SEARCH:
void NOMAD::Parameters::set_LH_SEARCH ( int p0 , int pi )
{
  _to_be_checked = true;
  _LH_search_p0  = (p0 <= 0 ) ? 0 : p0;
  _LH_search_pi  = (pi <= 0 ) ? 0 : pi;
}

// set_DIRECTION_TYPE (1/2):
void NOMAD::Parameters::set_DIRECTION_TYPE ( NOMAD::direction_type dt )
{
  _to_be_checked = true;
  if ( dt == NOMAD::UNDEFINED_DIRECTION || dt == NOMAD::NO_DIRECTION )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: DIRECTION_TYPE" );
  _direction_types.insert ( dt );
}

// set_DIRECTION_TYPE (2/2):
void NOMAD::Parameters::set_DIRECTION_TYPE ( const std::set<NOMAD::direction_type> & dt )
{
  std::set<NOMAD::direction_type>::const_iterator it , end = dt.end();
  for ( it = dt.begin() ; it != end ; ++it )
    set_DIRECTION_TYPE ( *it );
}

// set_SEC_POLL_DIR_TYPE (1/2):
void NOMAD::Parameters::set_SEC_POLL_DIR_TYPE ( NOMAD::direction_type dt )
{
  _to_be_checked = true;
  if ( dt == NOMAD::UNDEFINED_DIRECTION )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: SEC_POLL_DIR_TYPE" );
  _sec_poll_dir_types.insert ( dt );
}

// set_SEC_POLL_DIR_TYPE (2/2):
void NOMAD::Parameters::set_SEC_POLL_DIR_TYPE
( const std::set<NOMAD::direction_type> & dt ) {
  std::set<NOMAD::direction_type>::const_iterator it , end = dt.end();
  for ( it = dt.begin() ; it != end ; ++it )
    set_SEC_POLL_DIR_TYPE ( *it );
}

// set_RHO:
void NOMAD::Parameters::set_RHO ( const NOMAD::Double & rho )
{
  if ( !rho.is_defined() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , "invalid parameter: RHO" );
  _to_be_checked = true;
  _rho           = rho;
}

// set_H_MIN:
void NOMAD::Parameters::set_H_MIN ( const NOMAD::Double & h_min )
{
  if ( !h_min.is_defined() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , "invalid parameter: H_MIN" );
  _to_be_checked = true;
  _h_min         = h_min;
}

// set_H_MAX_0:
void NOMAD::Parameters::set_H_MAX_0 ( const NOMAD::Double & h_max )
{
  _to_be_checked = true;
  _h_max_0       = ( h_max.is_defined() ) ? h_max : NOMAD::INF;
}

// set_H_NORM:
void NOMAD::Parameters::set_H_NORM ( NOMAD::hnorm_type h_norm )
{
  _to_be_checked = true;
  _h_norm        = h_norm;
}

// set_SCALING (1/2):
void NOMAD::Parameters::set_SCALING ( int index , const NOMAD::Double & value )
{
  _to_be_checked = true;
  if ( index < 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: SCALING" );
  if ( index >= _scaling.size() )
    _scaling.resize ( index + 1 );
  
  _scaling[index] = value;
}

// set_SCALING (2/2):
void NOMAD::Parameters::set_SCALING ( const NOMAD::Point & s )
{
  _to_be_checked = true;
  _scaling       = s;
}

// set_FIXED_VARIABLE (1/3):
void NOMAD::Parameters::set_FIXED_VARIABLE ( int index , const NOMAD::Double & value )
{
  _to_be_checked = true;
  if ( index < 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: FIXED_VARIABLE" );
  if ( index >= _fixed_variables.size() )
    _fixed_variables.resize ( index + 1 );
  
  _fixed_variables[index] = value;
}

// set_FIXED_VARIABLE (2/3):
void NOMAD::Parameters::set_FIXED_VARIABLE ( int index )
{
  _to_be_checked = true;
  if ( index < 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: FIXED_VARIABLE (index < 0)" );
  if ( _x0s.empty() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: FIXED_VARIABLE (no starting point defined)" );

  if ( index >= _x0s[0]->size() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: FIXED_VARIABLE (incompatible starting point)" );

  if ( index >= _fixed_variables.size() )
    _fixed_variables.resize ( index + 1 );
  
  _fixed_variables[index] = (*_x0s[0])[index];
}

// set_FIXED_VARIABLE (2/3):
void NOMAD::Parameters::set_FIXED_VARIABLE ( const NOMAD::Point & fv )
{
  _to_be_checked   = true;
  _fixed_variables = fv;
}

// set_LOWER_BOUND (1/2):
void NOMAD::Parameters::set_LOWER_BOUND ( int index, const NOMAD::Double & value )
{
  _to_be_checked = true;
  if (index < 0)
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: LOWER_BOUND" );
  if ( index >= _lb.size() )
    _lb.resize(index+1);
  if ( !_lb[index].is_defined() || value > _lb[index] )
    _lb[index] = value;
}

// set_LOWER_BOUND (2/2):
void NOMAD::Parameters::set_LOWER_BOUND ( const NOMAD::Point & lb )
{
  _to_be_checked = true;
  _lb            = lb;
}

// set_UPPER_BOUND (1/2):
void NOMAD::Parameters::set_UPPER_BOUND ( int index, const NOMAD::Double & value )
{
  _to_be_checked = true;
  if (index < 0)
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: UPPER_BOUND" );
  if ( index >= _ub.size() )
    _ub.resize (index + 1);
  if ( !_ub[index].is_defined() || value < _ub[index] )
    _ub[index] = value;
}

// set_UPPER_BOUND (2/2):
void NOMAD::Parameters::set_UPPER_BOUND ( const NOMAD::Point & ub )
{
  _to_be_checked = true;
  _ub            = ub;
}

// set_SGTE_COST:
void NOMAD::Parameters::set_SGTE_COST ( int c )
{
  _to_be_checked = true;
  _sgte_cost = ( c > 0 ) ? c : -1;
}

// set_SGTE_EVAL_SORT:
void NOMAD::Parameters::set_SGTE_EVAL_SORT ( bool ses )
{
  _to_be_checked  = true;
  _sgte_eval_sort = ses;
}

// set_HAS_SGTET:
void NOMAD::Parameters::set_HAS_SGTE ( bool hs )
{
  _to_be_checked  = true;
  _has_sgte       = hs;
}

// set_OPT_ONLY_SGTE:
void NOMAD::Parameters::set_OPT_ONLY_SGTE ( bool oos )
{
  _to_be_checked  = true;
  _opt_only_sgte  = oos;
}

// set_SGTE_EXE:
void NOMAD::Parameters::set_SGTE_EXE ( const std::string & bb_exe   ,
				       const std::string & sgte_exe   )
{
  _to_be_checked     = true;
  _sgte_exe[bb_exe]  = sgte_exe;
}

// set_BB_EXE (1/3):
void NOMAD::Parameters::set_BB_EXE ( const std::string & bbexe )
{
  _to_be_checked = true;
  if ( _bb_output_type.empty() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: BB_EXE - BB_OUTPUT_TYPE must be defined first" );
  _bb_exe.clear();
  size_t nk = _bb_output_type.size();
  for ( size_t k = 0 ; k < nk ; ++k )
    _bb_exe.push_back ( bbexe );
}

// set_BB_EXE (2/3):
void NOMAD::Parameters::set_BB_EXE ( int m , const std::string * bbexe )
{
  _to_be_checked = true;

  if ( m <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: BB_EXE" );

  if ( m != static_cast<int>(_bb_output_type.size()) )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: BB_EXE - number of names or BB_OUTPUT_TYPE undefined" );

  size_t nk = _bb_output_type.size();
  for ( size_t k = 0 ; k < nk ; ++k )
    _bb_exe.push_back ( bbexe[k] );
}

// set_BB_EXE (3/3):
void NOMAD::Parameters::set_BB_EXE ( const std::list<std::string> & bbexe )
{
  _to_be_checked = true;
  if ( !bbexe.empty() && bbexe.size() != _bb_output_type.size() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: BB_EXE - number of names or BB_OUTPUT_TYPE undefined" );
  _bb_exe = bbexe;
}

// set_BB_INPUT_INCLUDE_TAG:
void NOMAD::Parameters::set_BB_INPUT_INCLUDE_TAG ( bool bbiit )
{
  _to_be_checked        = true;
  _bb_input_include_tag = bbiit;
}

// set_BB_INPUT_INCLUDE_SEED:
void NOMAD::Parameters::set_BB_INPUT_INCLUDE_SEED ( bool bbiis )
{
  _to_be_checked         = true;
  _bb_input_include_seed = bbiis;
}

// set_BB_REDIRECTION:
void NOMAD::Parameters::set_BB_REDIRECTION ( bool bbr )
{
  _to_be_checked  = true;
  _bb_redirection = bbr;
}

// set_BB_INPUT_TYPE (1/3):
void  NOMAD::Parameters::set_BB_INPUT_TYPE ( int index , NOMAD::bb_input_type bbit )
{
  _to_be_checked  = true;
  if ( index < 0 || index >= _dimension ||
       static_cast<int>(_bb_input_type.size()) != _dimension )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: BB_INPUT_TYPE" );
  _bb_input_type[index] = bbit;
}

// set_BB_INPUT_TYPE (2/3):
void NOMAD::Parameters::set_BB_INPUT_TYPE
( const std::vector<NOMAD::bb_input_type > & bbit )
{
  int n = bbit.size();
  for ( int i = 0 ; i < n ; ++i )
    set_BB_INPUT_TYPE ( i , bbit[i] );
}

// set_BB_INPUT_TYPE (3/3):
void NOMAD::Parameters::set_BB_INPUT_TYPE
( const std::list<NOMAD::bb_input_type > & bbit )
{
  int i = 0;
  std::list<NOMAD::bb_input_type>::const_iterator it , end = bbit.end();
  for ( it = bbit.begin() ; it != end ; ++it , ++i )
    set_BB_INPUT_TYPE ( i , *it );
}

// reset_PEB_changes:
void NOMAD::Parameters::reset_PEB_changes ( void ) const
{
  size_t nk = _bb_output_type.size();
  for ( size_t k = 0 ; k < nk ; ++k )
    if ( _bb_output_type[k] == NOMAD::PEB_E )
      _bb_output_type[k] = NOMAD::PEB_P;
}

// change_PEB_constraint_status:
void NOMAD::Parameters::change_PEB_constraint_status ( int index ) const
{
  if ( index < 0                                         ||
       index >= static_cast<int>(_bb_output_type.size()) ||
       _bb_output_type[index] != NOMAD::PEB_P               )
    throw NOMAD::Exception ( "Parameters.cpp" , __LINE__ ,
	  "error in Parameters::change_PEB_constraint_status(i): bad i" );
  _bb_output_type[index] = NOMAD::PEB_E;
}

// set_BB_OUTPUT_TYPE (1/2):
void NOMAD::Parameters::set_BB_OUTPUT_TYPE
( const std::list<NOMAD::bb_output_type> & bbot )
{
  int i = 0;
  std::vector<NOMAD::bb_output_type> bbot_vector ( bbot.size() );
  std::list<NOMAD::bb_output_type>::const_iterator end = bbot.end() , it;
  for ( it = bbot.begin() ; it != end ; ++it )
    bbot_vector[i++] = *it;
  set_BB_OUTPUT_TYPE ( bbot_vector );
}

// set_BB_OUTPUT_TYPE (2/2):
void NOMAD::Parameters::set_BB_OUTPUT_TYPE
( const std::vector<NOMAD::bb_output_type> & bbot )
{
  _to_be_checked          = true;

  _barrier_type           = NOMAD::EB;
  _has_constraints        = false;
  _has_EB_constraints     = false;
  _has_filter_constraints = false;

  _bb_output_type.clear();

  int m = static_cast<int>(bbot.size());

  if ( m <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: BB_OUTPUT_TYPE" );
  if ( !_bb_output_type.empty() &&
       m != static_cast<int>(_bb_output_type.size()) )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: BB_OUTPUT_TYPE - number of types" );

  _bb_output_type.resize (m);

  bool filter_used = false;
  bool pb_used     = false;
  bool peb_used    = false;

  _index_obj.clear();

  for ( int i = 0 ; i < m ; ++i ) {

    _bb_output_type[i] = bbot[i];

    switch ( bbot[i] ) {

    case NOMAD::OBJ:
      _index_obj.push_back(i);
      break;

    case NOMAD::EB:
      _has_constraints    = true;
      _has_EB_constraints = true;
      break;

    case NOMAD::FILTER:
      _has_constraints        = true;
      _has_filter_constraints = true;
      filter_used             = true;
      break;

    case NOMAD::PB:
      _has_constraints        = true;
      _has_filter_constraints = true;
      pb_used                 = true;
      break;
      
    case NOMAD::PEB_P:
    case NOMAD::PEB_E:
      _has_constraints        = true;
      _has_filter_constraints = true;
      pb_used                 = true;
      peb_used                = true;
      _bb_output_type[i]      = NOMAD::PEB_P;
      break;
    default:
      break;
    }
  }

  if ( _index_obj.empty() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: BB_OUTPUT_TYPE - OBJ not given" );
  if ( filter_used && pb_used )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: BB_OUTPUT_TYPE - F and PB/PEB used together" );

  if ( filter_used )
    _barrier_type = NOMAD::FILTER;
  else if ( pb_used )
    _barrier_type = (peb_used) ? NOMAD::PEB_P : NOMAD::PB;
}

// set_PROBLEM_DIR:
void NOMAD::Parameters::set_PROBLEM_DIR ( const std::string & dir )
{
  _to_be_checked = true;
  _problem_dir   = dir;
  if ( !_problem_dir.empty() && !NOMAD::Parameters::check_directory ( _problem_dir ) )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: PROBLEM_DIR" );
}

// set_TMP_DIR:
void NOMAD::Parameters::set_TMP_DIR ( const std::string & dir )
{
  _to_be_checked = true;
  _tmp_dir       = dir;
  if ( !_tmp_dir.empty() && !NOMAD::Parameters::check_directory ( _tmp_dir ) )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: TMP_DIR" );
}

// set_ADD_SEED_TO_FILE_NAMES:
void NOMAD::Parameters::set_ADD_SEED_TO_FILE_NAMES ( bool astfn )
{
  _to_be_checked          = true;
  _add_seed_to_file_names = astfn;
}

// set_SOLUTION_FILE:
void NOMAD::Parameters::set_SOLUTION_FILE ( const std::string & sf )
{
  _to_be_checked = true;
  _solution_file = sf;
  if ( sf.empty() )
    return;
  if ( !NOMAD::Parameters::check_directory ( _solution_file ) )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: SOLUTION_FILE" );
  _solution_file.resize ( _solution_file.size()-1 );
}

// set_HISTORY_FILE:
void NOMAD::Parameters::set_HISTORY_FILE ( const std::string & hf )
{
  _to_be_checked = true;
  _history_file  = hf;
  if ( hf.empty() )
    return;
  if ( !NOMAD::Parameters::check_directory ( _history_file ) )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: HISTORY_FILE" );
  _history_file.resize ( _history_file.size()-1 );
}

// set_CACHE_FILE:
void NOMAD::Parameters::set_CACHE_FILE ( const std::string & cf )
{
  _to_be_checked = true;
  _cache_file    = cf;
  if ( cf.empty() )
    return;
  if ( !NOMAD::Parameters::check_directory ( _cache_file ) )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: CACHE_FILE" );
  _cache_file.resize ( _cache_file.size()-1 );
}

// set_SGTE_CACHE_FILE:
void NOMAD::Parameters::set_SGTE_CACHE_FILE ( const std::string & cf )
{
  _to_be_checked   = true;
  _sgte_cache_file = cf;
  if ( cf.empty() )
    return;
  if ( !NOMAD::Parameters::check_directory ( _sgte_cache_file ) ) {
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: SGTE_CACHE_FILE");
  }
  _sgte_cache_file.resize ( _sgte_cache_file.size()-1 );
}

// set_X0:
// add a new point in the list of starting points:
void NOMAD::Parameters::set_X0 ( const NOMAD::Point & x0 )
{
  _to_be_checked = true;
  _x0s.push_back ( new NOMAD::Point ( x0 ) );
}

// indicate a x0 file or a cache file containing starting points:
void NOMAD::Parameters::set_X0 ( const std::string & file_name )
{
  _to_be_checked = true;

  if ( _dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "Parameters::set_X0() has been used before setting DIMENSION" );

  NOMAD::Point  tmp_x0 ( _dimension );
  std::string   complete_file_name = _problem_dir + file_name;
  std::ifstream fin ( complete_file_name.c_str() );
	
  if ( fin.fail() ) {
    std::string err = "invalid parameter: X0 - could not open file \'"
                      + complete_file_name + "\'";
    fin.close();
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ , err );
  }

  bool flag = true;
  try {
    fin >> tmp_x0;
  }
  catch ( NOMAD::Point::Bad_Input & ) {
    flag = false;

    // we suppose that the file name corresponds to a cache file:
    _x0_cache_file = file_name;
  }

  while ( flag ) {

    set_X0 ( tmp_x0 );

    // other starting points in the file ?
    flag = true;
    try {
      fin >> tmp_x0;
    }
    catch ( NOMAD::Point::Bad_Input & ) {
      flag = false;
    }
  }
  
  fin.close();
}

// set_DISPLAY_STATS (1/2):
void NOMAD::Parameters::set_DISPLAY_STATS ( const std::list<std::string> & ls )
{
  _display_stats = ls;
}

// set_DISPLAY_STATS (2/2):
void NOMAD::Parameters::set_DISPLAY_STATS ( const std::string & stats )
{
  if ( stats.empty() ) {
    _display_stats.clear();
    return;
  }

  NOMAD::Parameter_Entry pe ( "DISPLAY_STATS " + stats , false );

  std::list<std::string>::const_iterator end = pe.get_values().end() ,
                                         it  = pe.get_values().begin();
  std::list<std::string>                 ls;

  while ( it != end ) {
    ls.push_back ( *it );
    ++it;
  }

  ls.resize ( ls.size()-1 );

  set_DISPLAY_STATS ( ls );
}

// set_STATS_FILE (1/2):
void NOMAD::Parameters::set_STATS_FILE ( const std::string            & file_name ,
					 const std::list<std::string> & ls          )
{
  if ( file_name.empty() ) {
    reset_stats_file();
    return;
  }

  _to_be_checked   = true;
  _stats_file      = ls;
  _stats_file_name = file_name;

  if ( !NOMAD::Parameters::check_directory ( _stats_file_name ) )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: STATS_FILE" );

  _stats_file_name.resize ( _stats_file_name.size()-1 );
}

// set_STATS_FILE (2/2):
void NOMAD::Parameters::set_STATS_FILE ( const std::string & file_name ,
					 const std::string & stats       )
{
  NOMAD::Parameter_Entry pe ( "STATS_FILE " + file_name + " " + stats , false );
  std::list<std::string>::const_iterator end = pe.get_values().end() ,
                                         it  = pe.get_values().begin();
  std::list<std::string>                 ls;
  ++it;
  while ( it != end ) {
    ls.push_back(*it);
    ++it;
  }

  ls.resize ( ls.size()-1 );

  set_STATS_FILE ( file_name , ls );
}

// set_DISPLAY_DEGREE (1/3) - (accepts also dd_type arguments):
bool NOMAD::Parameters::set_DISPLAY_DEGREE ( int dd )
{
#ifndef DEBUG
  return set_DISPLAY_DEGREE ( NOMAD::itos(dd) );
#endif
  return true;
}

// set_DISPLAY_DEGREE (2/3):
bool NOMAD::Parameters::set_DISPLAY_DEGREE ( const std::string & dd )
{
#ifndef DEBUG
  {
    std::string ddu = dd;
    NOMAD::toupper ( ddu );
    
    if ( ddu == "NO" || ddu == "NO_DISPLAY" ) {
      set_DISPLAY_DEGREE ( 0 , 0 , 0 , 0 );
      return true;
    }

    else if ( ddu == "NORMAL" || ddu == "NORMAL_DISPLAY" ) {
      set_DISPLAY_DEGREE ( 1 , 1 , 1 , 1 );
      return true;
    }

    else if ( ddu == "FULL" || ddu == "FULL_DISPLAY" ) {
      set_DISPLAY_DEGREE ( 2 , 2 , 2 , 2 );
      return true;
    }
  }

  if ( dd.size() == 1 ) {
    int i;
    if ( !NOMAD::atoi ( dd[0] , i ) )
      return false;
    _out.set_degrees ( NOMAD::Display::int_to_dd(i) );
    return true;
  }

  if ( dd.size() != 4 )
    return false;

  // 1. general display:
  int gdd;
  if ( !NOMAD::atoi ( dd[0] , gdd ) )
    return false;

  // 2. search display:
  int sdd;
  if ( !NOMAD::atoi ( dd[1] , sdd ) )
    return false;

  // 3. poll display:
  int pdd;
  if ( !NOMAD::atoi ( dd[2] , pdd ) )
    return false;

  // 4. iterative display:
  int idd;
  if ( !NOMAD::atoi ( dd[3] , idd ) )
    return false;

  set_DISPLAY_DEGREE ( gdd , sdd , pdd , idd );

#endif
  return true;
}

// set_DISPLAY_DEGREE (3/3):
void NOMAD::Parameters::set_DISPLAY_DEGREE ( int gen_dd    ,
					     int search_dd ,
					     int poll_dd   ,
					     int iter_dd     )
{
#ifndef DEBUG
  _out.set_degrees ( NOMAD::Display::int_to_dd ( gen_dd    ) ,
		     NOMAD::Display::int_to_dd ( search_dd ) ,
		     NOMAD::Display::int_to_dd ( poll_dd   ) ,
		     NOMAD::Display::int_to_dd ( iter_dd   )   );
#endif
}

// set_OPEN_BRACE:
void NOMAD::Parameters::set_OPEN_BRACE ( const std::string & ob )
{
  _to_be_checked = true;
  _out.set_open_brace ( ob );
}

 // set_CLOSED_BRACE:
void NOMAD::Parameters::set_CLOSED_BRACE ( const std::string & cb )
{
  _to_be_checked = true;
  _out.set_closed_brace ( cb );
}

// set_SEED:
void NOMAD::Parameters::set_SEED ( int t )
{
  _to_be_checked = true;
  _seed          = ( t < 0 ) ? NOMAD::get_pid() : t;
}

// set_MAX_EVAL:
void NOMAD::Parameters::set_MAX_EVAL ( int e )
{
  _to_be_checked = true;
  _max_eval      = ( e <= 0 ) ? -1 : e;
}

// set_MAX_BB_EVAL:
void NOMAD::Parameters::set_MAX_BB_EVAL ( int bbe )
{
  _to_be_checked   = true;
  _max_bbe_decided = true;
  _max_bb_eval     = ( bbe < 0 ) ? -1 : bbe;
}

// set_MAX_SIM_BB_EVAL:
void NOMAD::Parameters::set_MAX_SIM_BB_EVAL ( int bbe )
{
  _to_be_checked   = true;
  _max_sim_bb_eval = ( bbe <= 0 ) ? -1 : bbe;
}

// set_MAX_SGTE_EVAL:
void NOMAD::Parameters::set_MAX_SGTE_EVAL ( int bbe )
{
  _to_be_checked = true;
  _sgte_max_eval = ( bbe < 0 ) ? -1 : bbe;
}

// set_MAX_TIME:
void NOMAD::Parameters::set_MAX_TIME ( int t )
{
  _to_be_checked = true;
  _max_time      = ( t <= 0 ) ? -1 : t;
}

// set_MAX_ITERATIONS:
void NOMAD::Parameters::set_MAX_ITERATIONS ( int it )
{
  _to_be_checked  = true;
  _max_iterations = ( it < 0 ) ? -1 : it;
}

// set_MAX_CACHE_MEMORY::
void NOMAD::Parameters::set_MAX_CACHE_MEMORY ( float mcm )
{
  _to_be_checked    = true;
  _max_cache_memory = ( mcm < 0.0 ) ? -1 : mcm;
}

// set_CACHE_SAVE_PERIOD:
void NOMAD::Parameters::set_CACHE_SAVE_PERIOD ( int csp )
{
  _to_be_checked     = true;
  _cache_save_period = ( csp <= 0 ) ? -1 : csp;
}

// set_STOP_IF_FEASIBLE:
void NOMAD::Parameters::set_STOP_IF_FEASIBLE ( bool sif )
{
  _to_be_checked    = true;
  _stop_if_feasible = sif;
}

// set_F_TARGET (1/2):
void NOMAD::Parameters::set_F_TARGET ( const NOMAD::Double & f_target )
{
  _to_be_checked = true;
  _f_target      = NOMAD::Point ( 1 , f_target );
}

// set_F_TARGET (2/2):
void NOMAD::Parameters::set_F_TARGET ( const NOMAD::Point & f_target )
{
  _to_be_checked = true;
  _f_target      = f_target;
}

// set_STAT_SUM_TARGET:
void NOMAD::Parameters::set_STAT_SUM_TARGET ( const NOMAD::Double & sst )
{
  _to_be_checked   = true;
  _stat_sum_target = sst;
}

// set_L_CURVE_TARGET:
void NOMAD::Parameters::set_L_CURVE_TARGET ( const NOMAD::Double & lct )
{
  _to_be_checked  = true;
  _L_curve_target = lct;
}

// set_MESH_UPDATE_BASIS:
void NOMAD::Parameters::set_MESH_UPDATE_BASIS ( const NOMAD::Double & mub )
{
  if ( !mub.is_defined() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: MESH_UPDATE_BASIS" );
  _to_be_checked     = true;
  _mesh_update_basis = mub;
}

// set_INITIAL_MESH_INDEX:
void NOMAD::Parameters::set_INITIAL_MESH_INDEX ( int ell_0 )
{
  _to_be_checked = true;
  if ( ell_0 > NOMAD::L_LIMITS )
    ell_0 = NOMAD::L_LIMITS;
  else if ( ell_0 < - NOMAD::L_LIMITS )
    ell_0 = -NOMAD::L_LIMITS;
  _initial_mesh_index = ell_0;
}

// set_MESH_COARSENING_EXPONENT:
void NOMAD::Parameters::set_MESH_COARSENING_EXPONENT ( int mce )
{
  _to_be_checked = true;
  if ( mce < 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: MESH_COARSENING_EXPONENT");
  _mesh_coarsening_exponent = mce;
}

// set_MESH_REFINING_EXPONENT:
void NOMAD::Parameters::set_MESH_REFINING_EXPONENT ( int mre )
{
  _to_be_checked = true;
  if ( mre >= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: MESH_REFINING_EXPONENT");
  _mesh_refining_exponent = mre;
}

// set_MAX_MESH_INDEX:
void  NOMAD::Parameters::set_MAX_MESH_INDEX ( int lmax )
{
  _to_be_checked = true;
  if ( lmax > NOMAD::L_LIMITS )
    lmax = NOMAD::L_LIMITS;
  _max_mesh_index = lmax;
}

// set_INITIAL_MESH_SIZE (1/3):
void NOMAD::Parameters::set_INITIAL_MESH_SIZE ( int                   index    ,
						const NOMAD::Double & d        ,
						bool                  relative   )
{
  if ( index < 0 || index >= _initial_mesh_size.size() || !d.is_defined() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: INITIAL_MESH_SIZE" );
  _to_be_checked = true;

  if ( relative ) {

    if ( !_lb.is_defined()  || !_ub.is_defined() )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: INITIAL_MESH_SIZE - bounds not defined" );
    
    if ( !_lb[index].is_defined()  || !_ub[index].is_defined() ||
	 d <= 0.0 || d > 1.0 )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: INITIAL_MESH_SIZE - relative value" );

    NOMAD::Double d2 = d;
    d2 *= _ub[index] - _lb[index];
    _initial_mesh_size[index] = d2;
  }
  else
    _initial_mesh_size[index] = d;
}

// set_INITIAL_MESH_SIZE (2/3):
void NOMAD::Parameters::set_INITIAL_MESH_SIZE ( const NOMAD::Point & delta_m_0 ,
						bool                 relative    )
{
  _to_be_checked = true;
  if ( relative ) {
    int nd = delta_m_0.size();
    for ( int i = 0 ; i < nd ; ++i )
      set_INITIAL_MESH_SIZE ( i , delta_m_0[i] , true );
  }
  else
    _initial_mesh_size = delta_m_0;
}

// set_INITIAL_MESH_SIZE (3/3):
void NOMAD::Parameters::set_INITIAL_MESH_SIZE ( const NOMAD::Double & d , bool relative )
{
  if ( _dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: INITIAL_MESH_SIZE - undefined dimension" );
  _to_be_checked = true;
  if ( relative )
    for ( int i = 0 ; i < _dimension ; ++i )
      set_INITIAL_MESH_SIZE ( i , d , true );
  else
    _initial_mesh_size = NOMAD::Point ( _dimension , d );
}

// set_MIN_MESH_SIZE (1/3):
void NOMAD::Parameters::set_MIN_MESH_SIZE ( int                   index    ,
					    const NOMAD::Double & d        ,
					    bool                  relative   )
{
  if (_dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: MIN_MESH_SIZE - undefined dimension" );
  
  if ( !_min_mesh_size.is_defined() )
    _min_mesh_size = NOMAD::Point ( _dimension );

  if ( index < 0 || index >= _min_mesh_size.size() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: MIN_MESH_SIZE" );
  _to_be_checked = true;
  if ( relative ) {
    
    if ( !_lb.is_defined()  || !_ub.is_defined() )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: MIN_MESH_SIZE - bounds not defined" );

    if ( !_lb[index].is_defined()  || !_ub[index].is_defined() ||
 	 !d.is_defined() || d <= 0.0 || d > 1.0 )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: MIN_MESH_SIZE - relative value" );
    NOMAD::Double d2 = d;
    d2 *= _ub[index] - _lb[index];
    _min_mesh_size[index] = d2;  
  }
  else
    _min_mesh_size[index] = d;
}

// set_MIN_MESH_SIZE (2/3):
void NOMAD::Parameters::set_MIN_MESH_SIZE ( const NOMAD::Point & delta_p_min ,
					    bool                 relative      )
{
  _to_be_checked = true;
  if ( relative ) {
    int nd = delta_p_min.size();
    for ( int i = 0 ; i < nd ; ++i )
      set_MIN_MESH_SIZE ( i , delta_p_min[i] , true );
  }
  else
    _min_mesh_size = delta_p_min;
}

// set_MIN_MESH_SIZE (3/3):
void NOMAD::Parameters::set_MIN_MESH_SIZE ( const NOMAD::Double & d , bool relative )
{
  if ( _dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: MIN_MESH_SIZE - undefined dimension" );
  _to_be_checked = true;
  if ( relative )
    for ( int i = 0 ; i < _dimension ; ++i )
      set_MIN_MESH_SIZE ( i , d , true );
  else
    _min_mesh_size = NOMAD::Point ( _dimension , d );
}

// set_MIN_POLL_SIZE (1/3):
void NOMAD::Parameters::set_MIN_POLL_SIZE ( int                   index    ,
					    const NOMAD::Double & d        ,
					    bool                  relative   )
{
  if ( _dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: MIN_POLL_SIZE - undefined dimension" );
  
  if ( !_min_poll_size.is_defined() )
    _min_poll_size = NOMAD::Point ( _dimension );

  if ( index < 0 || index >= _min_poll_size.size() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: MIN_POLL_SIZE" );
  _to_be_checked = true;
  if ( relative ) {
    if ( !_lb[index].is_defined()  || !_ub[index].is_defined() ||
 	 !d.is_defined() || d <= 0.0 || d > 1.0 )
      throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	    "invalid parameter: MIN_POLL_SIZE - relative value" );
    _min_poll_size[index] = d * ( _ub[index] - _lb[index] );
  }
  else
    _min_poll_size[index] = d;
}

// set_MIN_POLL_SIZE (2/3):
void NOMAD::Parameters::set_MIN_POLL_SIZE ( const NOMAD::Point & delta_p_min ,
					    bool                 relative      )
{
  _to_be_checked = true;
  if ( relative ) {
    int nd = delta_p_min.size();
    for ( int i = 0 ; i < nd ; ++i )
      set_MIN_POLL_SIZE ( i , delta_p_min[i] , true );
  }
  else
    _min_poll_size = delta_p_min;
}

// set_MIN_POLL_SIZE (3/3):
void NOMAD::Parameters::set_MIN_POLL_SIZE ( const NOMAD::Double & d , bool relative )
{
  if ( _dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: MIN_POLL_SIZE - undefined dimension" );
  _to_be_checked = true;
  if ( relative )
    for ( int i = 0 ; i < _dimension ; ++i )
      set_MIN_POLL_SIZE ( i , d , true );
  else
    _min_poll_size = NOMAD::Point ( _dimension , d );
}

// set_EXTENDED_POLL_TRIGGER:
void NOMAD::Parameters::set_EXTENDED_POLL_TRIGGER ( const NOMAD::Double & ept      ,
						    bool                  relative   )
{
  _to_be_checked = true;

  if ( !ept.is_defined() )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: EXTENDED_POLL_TRIGGER (undefined)" );

  if ( ept <= 0.0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: EXTENDED_POLL_TRIGGER: must be strictly positive" );
    
  _extended_poll_trigger = ept;
  _relative_ept          = relative;
}

// set_EXTENDED_POLL_ENABLED:
void NOMAD::Parameters::set_EXTENDED_POLL_ENABLED ( bool epe )
{
  _to_be_checked         = true;
  _extended_poll_enabled = epe;
}

// set_USER_CALLS_ENABLED:
void NOMAD::Parameters::set_USER_CALLS_ENABLED ( bool uce )
{
  _to_be_checked      = true;
  _user_calls_enabled = uce;
}

// set_ASYNCHRONOUS:
void NOMAD::Parameters::set_ASYNCHRONOUS ( bool a )
{
  _to_be_checked = true;
  _asynchronous  = a;
}

// set_MULTI_NB_MADS_RUNS:
void NOMAD::Parameters::set_MULTI_NB_MADS_RUNS ( int i )
{ 
  if ( i == 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: MULTI_NB_MADS_RUNS - has been set to zero" );

  _to_be_checked      = true;
  _multi_nb_mads_runs = ( i < 0 ) ? -1 : i;
}

// set_MULTI_OVERALL_BB_EVAL:
void NOMAD::Parameters::set_MULTI_OVERALL_BB_EVAL ( int i )
{
  _to_be_checked         = true;
  _multi_overall_bb_eval = ( i < 0 ) ? -1 : i;
}  

// set_MULTI_USE_DELTA_CRIT:
void NOMAD::Parameters::set_MULTI_USE_DELTA_CRIT ( bool b )
{
  _to_be_checked        = true;
  _multi_use_delta_crit = b;
}

// set_MULTI_F_BOUNDS:
void NOMAD::Parameters::set_MULTI_F_BOUNDS ( const NOMAD::Point & p )
{
  _to_be_checked  = true;
  if ( p.size() != 4 || p[0] >= p[1] || p[2] >= p[3] )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
			      "invalid parameter: MULTI_F_BOUNDS" );
  _multi_f_bounds = p;
}

// set_MULTI_FORMULATION:
void NOMAD::Parameters::set_MULTI_FORMULATION ( NOMAD::multi_formulation_type mft )
{
  _to_be_checked     = true;
  _multi_formulation = mft;
}

// set_OPPORTUNISTIC_LH:
void NOMAD::Parameters::set_OPPORTUNISTIC_LH ( bool opp )
{
  _to_be_checked     = true;
  _opportunistic_LH  = opp;
  _opp_LH_is_defined = true;
}

// set_OPPORTUNISTIC_CACHE_SEARCH:
void NOMAD::Parameters::set_OPPORTUNISTIC_CACHE_SEARCH ( bool opp )
{
  _to_be_checked              = true;
  _opportunistic_cache_search = opp;
  _opp_CS_is_defined          = true;
}

// set_OPPORTUNISTIC_EVAL:
void NOMAD::Parameters::set_OPPORTUNISTIC_EVAL ( bool oe )
{
  _to_be_checked      = true;
  _opportunistic_eval = oe;
}

// set_OPPORTUNISTIC_MIN_NB_SUCCESS:
void NOMAD::Parameters::set_OPPORTUNISTIC_MIN_NB_SUCCESS ( int mns )
{
  _to_be_checked = true;
  if ( mns <= 0 )
    mns = -1;
  _opportunistic_min_nb_success = mns;  
}

// set_OPPORTUNISTIC_MIN_EVAL:
void NOMAD::Parameters::set_OPPORTUNISTIC_MIN_EVAL ( int me )
{
  _to_be_checked = true;
  if ( me <= 0 )
    me = -1;
  _opportunistic_min_eval = me;
}

// set_OPPORTUNISTIC_MIN_F_IMPRVMT:
void NOMAD::Parameters::set_OPPORTUNISTIC_MIN_F_IMPRVMT ( const NOMAD::Double & mfi )
{
  _to_be_checked = true;
  if ( !mfi.is_defined() || mfi <= 0.0 )
    _opportunistic_min_f_imprvmt.clear();
  else
    _opportunistic_min_f_imprvmt = mfi;
}

// set_OPPORTUNISTIC_LUCKY_EVAL:
void NOMAD::Parameters::set_OPPORTUNISTIC_LUCKY_EVAL ( bool le )
{
  _to_be_checked            = true;
  _opportunistic_lucky_eval = le;
}

// set_PERIODIC_VARIABLE (1/2):
void NOMAD::Parameters::set_PERIODIC_VARIABLE ( int index )
{
  if ( _dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: PERIODIC_VARIABLE - undefined dimension" );

  if ( index < 0 || index >= _dimension )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: PERIODIC_VARIABLE - bad variable index" );

  if ( _periodic_variables.empty() )
    for ( int i = 0 ; i < _dimension ; ++i )
      _periodic_variables.push_back ( false );

  _periodic_variables[index] = true;
  _to_be_checked             = true;
}

// set_PERIODIC_VARIABLE (2/2):
void NOMAD::Parameters::set_PERIODIC_VARIABLE ( const std::vector<bool> & pv )
{
  _to_be_checked      = true;
  _periodic_variables = pv;
}

// set_VARIABLE_GROUP (1/3):
void NOMAD::Parameters::set_VARIABLE_GROUP ( const std::set<int> & var_indexes ,
					     int                   halton_seed   )
{  
  if ( _dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: VARIABLE_GROUP - undefined dimension" );

  if ( _bb_input_type.empty() ||
       static_cast<int>(_bb_input_type.size()) != _dimension )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: VARIABLE_GROUP - undefined blackbox input types" );

  _to_be_checked = true;
  
  std::set<NOMAD::direction_type> empty;

  _user_var_groups.insert ( new NOMAD::Variable_Group ( var_indexes ,
							empty       ,
							empty       ,
							halton_seed ,
							_out          ) );
}

// set_VARIABLE_GROUP (2/3):
void NOMAD::Parameters::set_VARIABLE_GROUP
( const std::set<int>                   & var_indexes        ,
  const std::set<NOMAD::direction_type> & direction_types    ,
  const std::set<NOMAD::direction_type> & sec_poll_dir_types ,
  int                                     halton_seed          )
{
  if ( _dimension <= 0 )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: VARIABLE_GROUP - undefined dimension" );

  if ( _bb_input_type.empty() ||
       static_cast<int>(_bb_input_type.size()) != _dimension )
    throw Invalid_Parameter ( "Parameters.cpp" , __LINE__ ,
	  "invalid parameter: VARIABLE_GROUP - undefined blackbox input types" );

  _to_be_checked = true;

  std::set<NOMAD::direction_type> dt = direction_types;
  if ( dt.empty() )
    dt.insert ( NOMAD::ORTHO_2N );

  _user_var_groups.insert ( new NOMAD::Variable_Group ( var_indexes        ,
							dt                 ,
							sec_poll_dir_types ,
							halton_seed        ,
							_out                 ) );
}

// set_VARIABLE_GROUP (3/3):
void NOMAD::Parameters::set_VARIABLE_GROUP
( const std::list<NOMAD::Variable_Group*> & vg )
{
  std::list<NOMAD::Variable_Group*>::const_iterator it , end = vg.end();
  for ( it = vg.begin() ; it != end ; ++it )
    set_VARIABLE_GROUP ( (*it)->get_var_indexes       () ,
 			 (*it)->get_direction_types   () ,
 			 (*it)->get_sec_poll_dir_types() ,
 			 (*it)->get_halton_seed       ()   );
}

/*-----------------------------------------------*/
/*  help (display parameter descriptions) - 1/3  */
/*-----------------------------------------------*/
void NOMAD::Parameters::help ( const std::string & param_name ) const
{
  std::list<std::string> ls;
  ls.push_back ( param_name );
  help ( ls );
}

/*-----------------------------------------------*/
/*  help (display parameter descriptions) - 2/3  */
/*-----------------------------------------------*/
void NOMAD::Parameters::help ( int argc , char ** argv ) const
{
  std::list<std::string> ls;
  if ( argc <= 2 )
    ls.push_back ( "ALL" );
  else
    for ( int i = 2 ; i < argc ; ++i )
      ls.push_back ( argv[i] );
  help ( ls );
}

/*-----------------------------------------------*/
/*  help (display parameter descriptions) - 3/3  */
/*-----------------------------------------------*/
void NOMAD::Parameters::help ( const std::list<std::string> & pnames ) const
{
#ifdef USE_MPI
  if ( !NOMAD::Slave::is_master() )
    return;
#endif

  // initial display:
  _out << std::endl
       << NOMAD::open_block ( "NOMAD - help on parameters - see more details in "
			      + NOMAD::USER_GUIDE_FILE );

  std::list<std::string> param_names = pnames;
  NOMAD::toupper ( param_names );

  bool   chk         = false;
  bool   display_all = false;

  if ( NOMAD::string_find ( "ALL PARAMS PARAMETERS EVERYTHING NOMAD" , param_names ) )
    display_all = true;

  // ADD_SEED_TO_FILE_NAMES:
  // -----------------------
  if ( display_all || NOMAD::string_find ( "ADD_SEED_TO_FILE_NAMES OUTPUTS \
                                            ADVANCED FILES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "ADD_SEED_TO_FILE_NAMES (advanced)" ) << std::endl
	 << ". if \'yes\', the seed is added to the name of" << std::endl
	 << "    output files (HISTORY_FILE, SOLUTION_FILE," << std::endl
	 << "    and STATS_FILE)" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'yes\'" << std::endl
	 << ". example: ADD_SEED_TO_FILE_NAMES no" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // ASYNCHRONOUS:
  // -------------
  if ( display_all || NOMAD::string_find ( "ASYNCHRONOUS PARALLELISM MPI \
                                            ADVANCED PMADS" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "ASYNCHRONOUS (advanced)" ) << std::endl
	 << ". asynchronous strategy for the parallel version" << std::endl
	 << ". if \'yes\', there can be evaluations in progress" << std::endl
	 << "    after an iteration has ended" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'yes\'" << std::endl
	 << ". example: ASYNCHRONOUS no" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // BB_EXE:
  // -------
  if ( display_all || NOMAD::string_find ( "BB_EXE BASIC BLACK-BOXES BLACKBOXES \
                                            EXECUTABLE FILES BI-OBJECTIVES \
                                            BIOBJECTIVES MULTIOBJECTIVES \
                                            BIMADS BI-MADS \
                                            MULTI-OBJECTIVES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "BB_EXE (basic)" ) << std::endl
	 << ". blackbox executable names" << std::endl
	 << ". list of strings" << std::endl
	 << ". no default, required (except in library mode)" << std::endl
	 << ". several executables can be used" << std::endl
	 << ". one executable can give several outputs" << std::endl
	 << ". use \' or \", and \'$\', to specify names or" << std::endl
	 << "    commands with spaces" << std::endl
	 << ". when the \'$\' character is put in first" << std::endl
	 << "    position of a string, it is considered" << std::endl
	 << "    as global and no path will be added" << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "BB_EXE bb.exe" << std::endl
	 << "BB_EXE bb1.exe bb2.exe" << std::endl
	 << "BB_EXE \'$nice bb.exe\'" << std::endl
	 << "BB_EXE \'$python bb.py\'" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }
  
  // BB_INPUT_INCLUDE_SEED:
  // ----------------------
  if ( display_all || NOMAD::string_find ( "BB_INPUT_INCLUDE_SEED BLACK-BOXES \
                                            BLACKBOXES \
                                            ADVANCED FILES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "BB_INPUT_INCLUDE_SEED (advanced)" ) << std::endl
	 << ". if \'yes\', the seed (\'SEED\') of the current" << std::endl
	 << "    execution is put as the first entry in" << std::endl
	 << "    all blackbox input files" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'no\'" << std::endl
	 << ". example: BB_INPUT_INCLUDE_SEED yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // BB_INPUT_INCLUDE_TAG:
  // ---------------------
  if ( display_all || NOMAD::string_find ( "BB_INPUT_INCLUDE_TAG \
                                            BLACK-BOXES BLACKBOXES \
                                            ADVANCED FILES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "BB_INPUT_INCLUDE_TAG (advanced)" ) << std::endl
	 << ". if \'yes\', the tag of a point is put as the first" << std::endl
	 << "    entry in all blackbox input files (second" << std::endl
	 << "    entry if BB_INPUT_INCLUDE_SEED=\'yes\')" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'no\'" << std::endl
	 << ". example: BB_INPUT_INCLUDE_TAG yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // BB_INPUT_TYPE:
  // --------------
  if ( display_all || NOMAD::string_find ( "BB_INPUT_TYPE BLACK-BOXES BLACKBOXES \
                                            BASIC VARIABLES \
                                            CONTINUOUS INTEGER BINARY" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "BB_INPUT_TYPE (basic)" ) << std::endl
	 << ". blackbox input types" << std::endl
	 << ". list of types for each variable" << std::endl
	 << ". " << NOMAD::open_block ( "available types" )
	 << "B: binary" << std::endl
	 << "C: categorical" << std::endl
	 << "I: integer" << std::endl
	 << "R: continuous" << std::endl
	 << NOMAD::close_block()
	 << ". default: * R (all continuous)" << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "BB_INPUT_TYPE ( R I B ) # for all 3 variables" << std::endl
	 << "BB_INPUT_TYPE 1-3 B     # variables 1 to 3 are binary" << std::endl
	 << "BB_INPUT_TYPE 0 I       # first variable is integer" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }

  // BB_OUTPUT_TYPE:
  // ---------------
  if ( display_all || NOMAD::string_find ( "BB_OUTPUT_TYPE BASIC CONSTRAINTS \
                                            BLACK-BOXES BLACKBOXES \
                                            BARRIER STATS BI-OBJECTIVES  h(x) \
                                            PB PEB FILTER STATS SUM AVG CNT_EVAL \
                                            COUNT PROGRESSIVE-BARRIER \
                                            MULTI-OBJECTIVE MINIMIZE \
                                            BIOBJECTIVES MULTIOBJECTIVES \
                                            BI-MADS BIMADS \
                                            OPTIMIZE" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "BB_OUTPUT_TYPE (basic)" ) << std::endl
	 << ". blackbox output types" << std::endl
	 << ". list of types for each blackbox output" << std::endl
	 << ". " << NOMAD::open_block ( "available types" ) << std::endl
	 << "OBJ       : objective value to minimize" << std::endl
	 << "            (define twice for bi-objective)" << std::endl
	 << "CSTR or PB: constraint <= 0 treated with" << std::endl
	 << "            Progressive Barrier (PB)" << std::endl
	 << "EB        : constraint <= 0 treated with" << std::endl
	 << "            Extreme Barrier (EB)" << std::endl
	 << "PEB       : constraint <= 0 treated with" << std::endl
	 << "            hybrid PB/EB" << std::endl
	 << "F         : constraint <= 0 treated with Filter" << std::endl
	 << "CNT_EVAL  : 0 or 1 output: count or not the" << std::endl
	 << "            evaluation" << std::endl
	 << "STAT_AVG  : NOMAD will compute the average" << std::endl
	 << "            value for this output"  << std::endl
	 << "STAT_SUM  : the same for the sum" << std::endl
	 << "NOTHING   : the output is ignored" << std::endl
	 << "-         : same as \'NOTHING\'" << std::endl
	 << NOMAD::close_block()
	 << ". STAT_SUM and STAT_AVG outputs have to be unique;" << std::endl
	 << "    they are updated at every new blackbox evaluation" << std::endl
	 << ". no default, required" << std::endl
	 << ". see user guide for blackbox output formats" << std::endl
	 << ". equality constraints are not natively supported" << std::endl
	 << ". see parameters LOWER_BOUND and UPPER_BOUND for" << std::endl
	 << "    bound constraints"    << std::endl
	 << ". " << NOMAD::open_block ( "examples" ) << std::endl
	 << "BB_EXE bb.exe                   # these two lines define" << std::endl
	 << "BB_OUTPUT_TYPE OBJ EB EB        # that bb.exe outputs" << std::endl
	 << "                                # three values" << std::endl << std::endl
	 << "BB_EXE bb1.exe bb2.exe bb2.exe  # bb1.exe outputs the" << std::endl
	 << "BB_OUTPUT_TYPE OBJ EB EB        # objective and bb2.exe" << std::endl
	 << "                                # two constraints" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }
  
  // BB_REDIRECTION:
  // ---------------
  if ( display_all || NOMAD::string_find ( "BB_REDIRECTION BLACK-BOXES BLACKBOXES \
                                            ADVANCED OUTPUT FILES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "BB_REDIRECTION (advanced)" ) << std::endl
	 << ". if NOMAD uses a redirection (\'>\') to" << std::endl
	 << "    create blackbox output files" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'yes\'" << std::endl
	 << ". if \'no\', the blackbox has to manage its" << std::endl
	 << "    own output files (see user guide)" 	 << std::endl
	 << ". example: BB_INPUT_INCLUDE_TAG yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // CACHE_FILE:
  // -----------
  if ( display_all || NOMAD::string_find ( "CACHE_FILE BASIC \
                                            FILES OUTPUTS" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "CACHE_FILE (basic)" ) << std::endl
	 << ". cache file; if the specified file" << std::endl
	 << "    does not exist, it will be created" << std::endl
	 << ". argument: one string" << std::endl
	 << ". no default" << std::endl
	 << ". points already in the file will be" << std::endl
	 << "    tagged as true evaluations" << std::endl
	 << ". example: CACHE_FILE cache.bin" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // CACHE_SAVE_PERIOD:
  // ------------------
  if ( display_all || NOMAD::string_find ( "CACHE_SAVE_PERIOD OUTPUTS \
                                            ITERATIONS ADVANCED FILES" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "CACHE_SAVE_PERIOD (advanced)" ) << std::endl
	   << ". period (iterations) at which the cache file is saved" << std::endl
	   << "    (if CACHE_FILE is defined; disabled for bi-objective)" << std::endl
	   << ". argument: one nonnegative integer" << std::endl
	   << ". default: 25" << std::endl
	   << ". example: CACHE_SAVE_PERIOD 10" << std::endl
	   << NOMAD::close_block();
      chk = true;
    }

  // CACHE_SEARCH:
  // -------------
  if ( display_all || NOMAD::string_find ( "CACHE_SEARCH ADVANCED\
                                            CACHE_FILE SGTE_CACHE_FILE" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "CACHE_SEARCH (advanced)" ) << std::endl
	   << ". enable or disable the cache search" << std::endl
	   << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	   << ". default: \'no\'" << std::endl
	   << ". the search looks in the cache between iterations" << std::endl
	   << ". this can be useful when a non-empty initial cache" << std::endl
	   << "    file is provided or with an extern cache that" << std::endl
	   << "    the user updates independently" << std::endl
	   << ". example: CACHE_SEARCH yes" << std::endl
	   << NOMAD::close_block();
      chk = true;
    }
  
  // CLOSED_BRACE:
  // -------------
  if ( display_all || NOMAD::string_find ( "CLOSED_BRACES INDENTATION TABULATIONS \
                                            BLOCKS ADVANCED DISPLAY" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "CLOSED_BRACE (advanced)" ) << std::endl
	   << ". string displayed at the end of indented" << std::endl
	   << "    blocks in full display mode" << std::endl
	   << ". argument: one string" << std::endl
	   << ". default: \'}\'" << std::endl
	   << ". example: CLOSED_BRACE End" << std::endl
	   << NOMAD::close_block();
      chk = true;
    }

  // DIMENSION:
  // ----------
  if ( display_all || NOMAD::string_find ( "DIMENSION BASIC VARIABLE" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "DIMENSION (basic)" ) << std::endl
	   << ". number of variables" << std::endl
	   << ". argument: one positive integer" << std::endl
	   << ". no default, required" << std::endl
	   << ". example: DIMENSION 3" << std::endl
	   << NOMAD::close_block();
      chk = true;
    }

  // DIRECTION_TYPE:
  // ---------------
  if ( display_all || NOMAD::string_find ( "DIRECTION_TYPES DIRECTIONS_TYPES \
                                            BASIC ORTHO LT GPS \
                                            MADS DIRECTIONS 2N N+1 POLL \
                                            ORTHOMADS ORTHO-MADS LTMADS LT-MADS \
                                            RANDOM STATIC UNIFORM ANGLES" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "DIRECTION_TYPE (basic)" ) << std::endl
	 << ". types of directions used in the poll step" << std::endl
	 << ". arguments: direction types (see user" << std::endl
	 << "             guide for available types)" << std::endl
	 << ". default: ORTHO (2n OrthoMADS directions)" << std::endl
	 << ". several direction types can be defined" << std::endl
	 << "    at the same time (one direction per line)" << std::endl
	 << ". " << NOMAD::open_block ( "examples" ) << std::endl
	 << "DIRECTION_TYPE ORTHO 1                  #  OrthoMADS, 1" << std::endl
	 << "DIRECTION_TYPE ORTHO 2                  #  OrthoMADS, 2" << std::endl
	 << "DIRECTION_TYPE ORTHO                    #  OrthoMADS, 2n" << std::endl
	 << "DIRECTION_TYPE ORTHO 2N                 #  OrthoMADS, 2n" << std::endl
	 << "DIRECTION_TYPE LT    1                  #  LT-MADS, 1" << std::endl
	 << "DIRECTION_TYPE LT    2                  #  LT-MADS, 2" << std::endl
	 << "DIRECTION_TYPE LT    N+1                #  LT-MADS, n+1" << std::endl
	 << "DIRECTION_TYPE LT                       #  LT-MADS, 2n" << std::endl
	 << "DIRECTION_TYPE LT    2N                 #  LT-MADS, 2n" << std::endl
	 << "DIRECTION_TYPE GPS   BINARY or BIN      #  GPS, bin var" << std::endl
	 << "DIRECTION_TYPE GPS   N+1                #  GPS, n+1," << std::endl
	 << "                                        #  static" << std::endl
	 << "DIRECTION_TYPE GPS   N+1 STATIC         #  GPS, n+1," << std::endl
	 << "                                        #  static" << std::endl
	 << "DIRECTION_TYPE GPS   N+1 STATIC UNIFORM #  GPS, n+1," << std::endl
	 << "                                        #  static, unif" << std::endl
	 << "DIRECTION_TYPE GPS   N+1 RAND           #  GPS, n+1," << std::endl
	 << "                                        #  rand" << std::endl
	 << "DIRECTION_TYPE GPS   N+1 RAND   UNIFORM #  GPS, n+1," << std::endl
	 << "                                        #  rand, unif" << std::endl
	 << "DIRECTION_TYPE GPS                      #  GPS, 2n," << std::endl
	 << "                                        #  static" << std::endl
	 << "DIRECTION_TYPE GPS   2N                 #  GPS, 2n," << std::endl
	 << "                                        #  static" << std::endl
	 << "DIRECTION_TYPE GPS   2N  STATIC         #  GPS, 2n," << std::endl
	 << "                                        #  static" << std::endl
	 << "DIRECTION_TYPE GPS   2N  RAND           #  GPS, 2n," << std::endl
	 << "                                        #  rand" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }
  
  // DISPLAY_DEGREE:
  // ---------------
  if ( display_all || NOMAD::string_find ( "DISPLAY_DEGREES OUTPUTS \
                                            SCREEN DEBUG" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "DISPLAY_DEGREE (basic)" ) << std::endl
	 << ". " << NOMAD::open_block ( "display degree" )
	 << "0: no display" << std::endl
	 << "1: normal display" << std::endl
	 << "2: full display" << std::endl
      	 << NOMAD::close_block()
	 << ". argument: one integer in { 0, 1, 2 } (basic)" << std::endl
	 << "         or one string in { \'NO\', \'NO_DISPLAY\', \'NORMAL\',"
	 << std::endl
	 << "                            \'NORMAL_DISPLAY\', \'FULL\'," << std::endl
	 << "                            \'FULL_DISPLAY\' }" << std::endl
	 << "         or one string composed of 4 integers each in" << std::endl
	 << "            { 0, 1, 2 } (advanced)" << std::endl
	 << ". default: 1" << std::endl
	 << ". "
	 << NOMAD::open_block("advanced use with 4 digits (see user guide for details)")
	 << "#1 general display degree   (before and after" << std::endl
	 << "                             all iterations)" << std::endl
	 << "#2 search display degree    (during searches)" << std::endl
	 << "#3 poll display degree      (during polls)" << std::endl
	 << "#4 iterative display degree (other displays at" << std::endl
	 << "                             each iteration)" << std::endl
	 << NOMAD::close_block()
	 << ". " << NOMAD::open_block ( "examples" )
	 << "DISPLAY_DEGREE 1    # basic   : normal display" << std::endl
	 << "DISPLAY_DEGREE 0020 # advanced: display only" << std::endl
	 << "                    #           poll info" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }

  // DISPLAY_STATS:
  // ---------------
  if ( display_all || NOMAD::string_find ( "DISPLAY_STATS OUTPUTS LATEX FORMAT \
                                            SCREEN DEBUG" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "DISPLAY_STATS (basic)" ) << std::endl
	 << ". format of the outputs displayed at each success" << std::endl
	 << "  (single-objective)" << std::endl
	 << ". format of the final Pareto front" << std::endl
	 << "  (multi-objective)" << std::endl
	 << ". arguments: list of strings possibly including" << std::endl
	 << "    the following keywords:" << std::endl
	 << "      BBE       : blackbox evaluations" << std::endl
	 << "      BBO       : blackbox outputs" << std::endl
	 << "      EVAL      : evaluations (includes cache hits)" << std::endl
      	 << "      MESH_INDEX: mesh index" << std::endl
	 << "      MESH_SIZE : mesh size delta_k^m" << std::endl
	 << "      OBJ       : objective function value" << std::endl
	 << "      POLL_SIZE : poll size delta_k^p" << std::endl
	 << "      SIM_BBE   : simulated blackbox evaluations" << std::endl
	 << "      SGTE      : surrogate evaluations" << std::endl
	 << "      SOL       : current feasible iterate" << std::endl
	 << "      STAT_SUM  : value of stat SUM" << std::endl
	 << "      STAT_AVG  : value of stat AVG" << std::endl
	 << "      TIME      : real time in seconds" << std::endl
	 << "      VARi      : value of variable i" << std::endl
	 << "                  (0 for the first variable)"
	 << std::endl
	 << ". all outputs may be formatted using C style" << std::endl
	 << "  (%f, %e, %E, %g, %G, %i, %d with the possibility" << std::endl
	 << "  to specify the display width and the precision)" << std::endl
	 << "  example: %5.2Ef displays f in 5 columns and 2 decimals" << std::endl
	 << "           in scientific notation" << std::endl
	 << ". do not use quotes" << std::endl
	 << ". the \'%\' character may be explicitely indicated with \'\\%\'"
	 << std::endl
	 << ". see details in user guide" << std::endl
	 << ". defaults: BBE OBJ (single-objective)" << std::endl
	 << "            OBJ     (multi-objective)"  << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "DISPLAY_STATS TIME f=OBJ" << std::endl
	 << "DISPLAY_STATS %5.2obj" << std::endl
	 << "# for LaTeX tables:" << std::endl
	 << "DISPLAY_STATS $BBE$ & ( $%12.5SOL, ) & $OBJ$ \\\\" << std::endl
      	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }

  // EPSILON:
  // --------
  if ( display_all || NOMAD::string_find ( "EPSILON ADVANCED \
                                            PRECISION REALS COMPARISONS" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "EPSILON (advanced)" ) << std::endl
	 << ". precision on reals" << std::endl
	 << ". argument: one positive real" << std::endl
	 << ". default: 1E-13" << std::endl
	 << ". example: EPSILON 1E-8" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // EXTENDED_POLL_ENABLED:
  // ----------------------
  if ( display_all || NOMAD::string_find ( "EXTENDED_POLL_ENABLED \
                                            EXTENDED_POLL_DISABLED \
                                            ADVANCED" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "EXTENDED_POLL_ENABLED (advanced)" ) << std::endl
	 << ". if \'no\', the extended poll for categorical" << std::endl
	 << "    variables is disabled" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'yes\'" << std::endl
	 << ". the extended poll uses the surrogate" << std::endl
	 << "    if HAS_SGTE or SGTE_EXE is defined" << std::endl
	 << ". example: EXTENDED_POLL_ENABLED yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // EXTENDED_POLL_TRIGGER:
  // ----------------------
  if ( display_all || NOMAD::string_find ( "EXTENDED_POLL_TRIGGER ADVANCED \
                                            MIXED MVP CATEGORICAL" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "EXTENDED_POLL_TRIGGER (advanced)" ) << std::endl
	 << ". extended poll trigger for categorical variables" << std::endl
	 << ". argument: one positive real (can be relative)" << std::endl
	 << ". an extended poll around the extended poll point y" << std::endl
	 << "    constructed from an iterate xk is performed if" << std::endl
	 << "    f(y) < f(xk)+trigger or f(y) < f(xk)+|f(x_k)|*trigger" << std::endl
	 << "    (relative value)" << std::endl
	 << ". see details in user guide" << std::endl
	 << ". default: r0.1" << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "EXTENDED_POLL_TRIGGER 10.0  # ext poll trigger of 10" << std::endl
	 << "EXTENDED_POLL_TRIGGER r0.2  # ext poll trigger of 20%" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }

  // F_TARGET:
  // ---------
  if ( display_all || NOMAD::string_find ( "F_TARGET BASIC BLACK-BOXES BLACKBOXES \
                                            BI-OBJECTIVES STOPPING \
                                            MULTI-OBJECTIVES \
                                            BIOBJECTIVES MULTIOBJECTIVES \
                                            BI-MADS BIMADS \
                                            TERMINATES TERMINATION" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "F_TARGET (basic)" ) << std::endl
	 << ". NOMAD terminates if fi <= F_TARGET[i] for" << std::endl
	 << "    all objectives i" << std::endl
	 << ". arguments: one or two reals (single or bi-obj.)" << std::endl
	 << ". no default" << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "F_TARGET 0.0         # single-objective" << std::endl
	 << "F_TARGET ( 0.0 0.0 ) # bi-objective" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }
  
  // FIXED_VARIABLE:
  // ---------------
  if ( display_all || NOMAD::string_find ( "FIXED_VARIABLE ADVANCED \
                                            FILES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "FIXED_VARIABLE (advanced)" ) << std::endl
	 << ". fix some variables to some specific values" << std::endl
	 << ". arguments: variable indexes and values" << std::endl
	 << ". no default" << std::endl
	 << ". values are optional if at least one starting point" << std::endl
	 << "    is defined" << std::endl
	 << ". can be given by a text file containing DIMENSION" << std::endl
	 << "    entrie (use \'-\' for free variables)" << std::endl
	 << ". variables inside groups defined by VARIABLE_GROUP" << std::endl
	 << "    cannot be fixed" << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "FIXED_VARIABLE ( 0.0 - 0.0 ) # variables 0 and 2 are fixed" << std::endl
	 << "FIXED_VARIABLE fixed.txt     # with a file" << std::endl
	 << "FIXED_VARIABLE 0-1 3.14      # 2 first variables fixed" << std::endl
	 << "                             # to 3.14" << std::endl
	 << "FIXED_VARIABLE 0             # first variable fixed to" << std::endl
	 << "                             # its X0 value" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }

  // H_MAX_0:
  // --------
  if ( display_all || NOMAD::string_find ( "H_MAX_0 HMAX_0 HMAX ADVANCED \
                                            CONSTRAINTS PB FILTER PEB \
                                            L1 L2 LINF L_INF \
                                            FEASIBILITY \
                                            PROGRESSIVE-BARRIER" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "H_MAX_0 (advanced)" ) << std::endl
	 << ". initial value of h_max (for PB and" << std::endl
	 << "    F constraints handling strategies)" << std::endl
	 << ". argument: one positive real" << std::endl
	 << ". default: 1E+20" << std::endl
	 << ". points x such that h(x) > h_max are" << std::endl
	 << "    rejected (h measures the feasibility)" << std::endl
	 << ". example: H_MAX_0 100.0" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // H_MIN:
  // ------
  if ( display_all || NOMAD::string_find ( "H_MIN HMIN ADVANCED \
                                            CONSTRAINTS PB FILTER PEB \
                                            L1 L2 LINF L_INF \
                                            FEASIBILITY \
                                            PROGRESSIVE-BARRIER" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "H_MIN (advanced)" ) << std::endl
	 << ". value of h_min; x is feasible if h(x) <= h_min" << std::endl
	 << "  (h measures the feasibility)" << std::endl
	 << ". argument: one positive real" << std::endl
	 << ". default: 0.0" << std::endl
	 << ". example: H_MIN 1E-5" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // H_NORM:
  // -------
  if ( display_all || NOMAD::string_find ( "H_NORM ADVANCED \
                                            CONSTRAINTS PB FILTER PEB \
                                            L1 L2 LINF L_INF \
                                            FEASIBILITY \
                                            PROGRESSIVE-BARRIER" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "H_NORM (advanced)" ) << std::endl
	 << ". norm used by the F and PB constraints handling" << std::endl
	 << "    strategies to compute h(x) (h measures the" << std::endl
	 << "    feasibility)" << std::endl
	 << ". argument: one string in { \'L1\' , \'L2\' , \'Linf\' }" << std::endl
	 << ". default: \'L2\'" << std::endl
	 << ". example: H_NORM Linf" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // HALTON_SEED:
  // ------------
  if ( display_all || NOMAD::string_find ( "HALTON_SEED BASIC ORTHO-MADS ORTHOMADS \
                                            RANDOM \
                                            DETERMINISTIC ORTHOGONAL DIRECTIONS \
                                            POLL" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "HALTON_SEED (basic)" ) << std::endl
	 << ". Halton seed to get different runs with OrthoMADS" << std::endl
	 << ". argument: one nonnegative integer" << std::endl
	 << ". default: the DIMENSION^th prime number" << std::endl
	 << ". choose values greater to or equal to this default" << std::endl
	 << ". see also parameter SEED for different LT-MADS runs" << std::endl
	 << ". example: HALTON_SEED 10" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // HAS_SGTE:
  // ---------
  if ( display_all || NOMAD::string_find ( "HAS_SGTE SGTE_EXE ADVANCED SURROGATES \
                                            BLACK-BOXES	BLACKBOXES \
                                            SGTES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "HAS_SGTE (advanced)" ) << std::endl
	 << ". to indicate that the problem has a surrogate" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'no\' if parameter SGTE_EXE is undefined," << std::endl
	 << "           \'yes\' otherwise" << std::endl
	 << ". this parameter is not necessary in batch" << std::endl
	 << "    mode, but essential in library mode when" << std::endl
	 << "    no surrogate executable is provided" << std::endl
	 << ". example: HAS_SGTE yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // HISTORY_FILE:
  // -------------
  if ( display_all || NOMAD::string_find ( "HISTORY_FILE BASIC \
                                            FILES OUTPUTS" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "HISTORY_FILE (basic)" ) << std::endl
	 << ". history file: contains all trial points" << std::endl
	 << ". includes multiple evaluations" << std::endl
	 << ". argument: one string" << std::endl
	 << ". no default" << std::endl
	 << ". the seed is added to the file name" << std::endl
	 << "    if ADD_SEED_TO_FILE_NAMES=\'yes\'" << std::endl
	 << ". example: HISTORY_FILE his.txt" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // INF_STR:
  // --------
  if ( display_all || NOMAD::string_find ( "INF_STR ADVANCED \
                                            INFINITY DISPLAY REALS" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "INF_STR (advanced)" ) << std::endl
	 << ". string used to display infinity" << std::endl
	 << ". argument: one string" << std::endl
	 << ". default: \"inf\"" << std::endl
	 << ". example: INF_STR Infinity" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // INITIAL_MESH_INDEX:
  // -------------------
  if ( display_all || NOMAD::string_find ( "INITIAL_MESH_INDEX ADVANCED \\DELTA \
                                            MADS L0 L_0 \\ELL_0" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "INITIAL_MESH_INDEX (advanced)" ) << std::endl
	 << ". initial mesh index \\ell_0" << std::endl
	 << ". argument: one integer (can be negative)" << std::endl
	 << ". default: 0" << std::endl
	 << ". example: INITIAL_MESH_INDEX -1" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // INITIAL_MESH_SIZE:
  // ------------------
  if ( display_all || NOMAD::string_find ( "INITIAL_MESH_SIZE BASIC MADS \
                                            \\DELTA_0 " , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "INITIAL_MESH_SIZE (basic)" ) << std::endl
	 << ". initial mesh size" << std::endl
	 << ". arguments: one or DIMENSION positive real(s)" << std::endl
	 << ". defaults: r0.1 if bounds are defined (10% of the range)," << std::endl
	 << "            1.0 otherwise" << std::endl
	 << ". NOMAD uses one mesh size per variable to achieve scaling" << std::endl
	 << ". values can be given with \'r\' to indicate a proportion of" << std::endl
	 << "    the bounds range (bounds have to be defined for the" << std::endl
	 << "    corresponding variables)" << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "INITIAL_MESH_SIZE 1.0          # for all variables" << std::endl
	 << "INITIAL_MESH_SIZE ( 3 - r0.1 ) # for all variables" << std::endl
	 << "                               # (default considered" << std::endl
	 << "                               # for 2nd variable)" << std::endl
	 << "INITIAL_MESH_SIZE 1 0.5        # for var. 1 only" << std::endl
	 << "INITIAL_MESH_SIZE 2-4 r0.25    # for var. 2 to 4" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }
  
  // L_CURVE_TARGET:
  // ---------------
  if ( display_all || NOMAD::string_find ( "L_CURVE_TARGET ADVANCED TERMINATION \
                                            STOPPING TERMINATES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "L_CURVE_TARGET (advanced)" ) << std::endl
	 << ". MADS terminates if it detects that the objective will" << std::endl
	 << "    not reach this value (based on an approximation" << std::endl
	 << "    of the L-shaped curve obj_value v.s. bb_eval)" << std::endl
	 << ". argument: one real" << std::endl
	 << ". no default" << std::endl
	 << ". example: L_CURVE_TARGET 10.0" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // LH_SEARCH:
  // ----------
  if ( display_all || NOMAD::string_find ( "LH_SEARCH LATIN-HYPERCUBE \
                                            SAMPLING BASIC" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "LH_SEARCH (basic)" ) << std::endl
	 << ". Latin-Hypercube sampling (search)" << std::endl
	 << ". arguments: two nonnegative integers p0 and pi" << std::endl
	 << ". defaults: no search for single-objective" << std::endl
	 << "         or one initial search for bi-objective" << std::endl
	 << "            (see user guide)" << std::endl
	 << ". p0: number of initial LH search points" << std::endl
	 << "      (or in first MADS run for bi-obj.)" << std::endl
	 << ". pi: LH search points at each iteration" << std::endl
	 << "      (or in 2nd MADS run for bi-obj.)" << std::endl
	 << ". the search can be opportunistic or not" << std::endl
	 << "    (see parameter OPPORTUNISTIC_LH)" << std::endl
	 << ". example: LH_SEARCH 100 0" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // LOWER_BOUND:
  // ------------
  if ( display_all || NOMAD::string_find ( "LOWER_BOUND BASIC LB BOUNDS FILES" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "LOWER_BOUND (basic)" ) << std::endl
	 << ". lower bounds for each variable" << std::endl
	 << ". no default" << std::endl
	 << ". can be defined by various methods (see user guide)" << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "LOWER_BOUND * 0.0   # all variables are nonnegative" << std::endl
	 << "LOWER_BOUND 0-2 0.0 # the 3 first var. are nonnegative" << std::endl
	 << "LOWER_BOUND 0 0.0   # the first var. is nonnegative" << std::endl
	 << "LOWER_BOUND lb.txt  # bounds are defined in \'lb.txt\'" << std::endl
	 << "                    # containing DIMENSION values" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }

  // MAX_BB_EVAL:
  // ------------
  if ( display_all || NOMAD::string_find ( "MAX_BB_EVAL BASIC BLACK-BOXES BLACKBOXES \
                                            MAXIMUM CACHE BBEVAL \
                                            NUMBER EVALUATIONS STOPPING \
                                            BI-0BJECTIVE MULTI-OBJECTIVE \
                                            BIOBJECTIVES MULTIOBJECTIVES \
                                            BI-MADS BIMADS \
                                            TERMINATES TERMINATION" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MAX_BB_EVAL (basic)" ) << std::endl
	 << ". maximum number of blackbox evaluations" << std::endl
	 << ". argument: one positive integer" << std::endl
	 << ". no default" << std::endl
	 << ". doesn\'t consider evaluations taken in the" << std::endl
	 << "    cache (cache hits)" << std::endl
	 << ". in bi-objective mode: max number of blackbox" << std::endl
	 << "    evaluations for each MADS run" << std::endl
	 << ". example: MAX_BB_EVAL 1000" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MAX_CACHE_MEMORY:
  // -----------------
  if ( display_all || NOMAD::string_find ( "MAX_CACHE_MEMORY ADVANCED \
                                            MAXIMUM RAM STOPPING \
                                            MB MEGA-BYTES MEGABYTES \
                                            TERMINATES TERMINATION" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MAX_CACHE_MEMORY (advanced)" ) << std::endl
	 << ". the program terminates as soon as the cache" << std::endl
	 << "    reaches this memory limit" << std::endl
	 << ". argument: one positive integer (expressed in MB)" << std::endl
	 << ". no default" << std::endl
	 << ". example: MAX_CACHE_MEMORY 1024 # limit of 1GB memory" << std::endl
	 << "                                 # occupation" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MAX_EVAL:
  // ---------
  if ( display_all || NOMAD::string_find ( "MAX_EVAL ADVANCED BLACK-BOXES BLACKBOXES \
                                            MAXIMUM CACHE BBEVAL \
                                            NUMBER EVALUATIONS STOPPING \
                                            TERMINATES TERMINATION" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MAX_EVAL (advanced)" ) << std::endl
	 << ". maximum number of evaluations" << std::endl
	 << ". argument: one positive integer" << std::endl
	 << ". no default" << std::endl
	 << ". includes evaluations taken in" << std::endl
	 << "    the cache (cache hits)" << std::endl
	 << ". example: MAX_EVAL 1000" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MAX_ITERATIONS:
  // ---------------
  if ( display_all || NOMAD::string_find ( "MAX_ITERATIONS ADVANCED \
                                            MAXIMUM MADS \
                                            NUMBER STOPPING \
                                            TERMINATES TERMINATION" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MAX_ITERATIONS (advanced)" ) << std::endl
	 << ". maximum number of MADS iterations" << std::endl
	 << ". argument: one positive integer" << std::endl
	 << ". no default" << std::endl
	 << ". example: MAX_ITERATIONS 20" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MAX_MESH_INDEX:
  // ---------------
  if ( display_all || NOMAD::string_find ( "MAX_MESH_INDEX ADVANCED MAXIMUM \
                                            \\DELTA TERMINATION \
                                            STOPPING TERMINATES \\ELL" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "MAX_MESH_INDEX (advanced)" ) << std::endl
	   << ". maximum mesh index" << std::endl
	   << ". argument: one nonnegative integer" << std::endl
	   << ". no default" << std::endl
	   << ". example: MAX_MESH_INDEX 5" << std::endl
	   << NOMAD::close_block();
      chk = true;
    }
  
  // MAX_SGTE_EVAL:
  // --------------
  if ( display_all || NOMAD::string_find ( "MAX_SGTE_EVAL ADVANCED BLACK-BOXES \
                                            BLACKBOXES \
                                            MAXIMUM SURROGATES BBEVAL SGTES \
                                            NUMBER EVALUATIONS STOPPING \
                                            TERMINATES TERMINATION" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MAX_SGTE_EVAL (advanced)" ) << std::endl
	 << ". maximum number of surrogate evaluations" << std::endl
	 << ". argument: one positive integer" << std::endl
	 << ". no default" << std::endl
	 << ". example: MAX_SGTE_EVAL 10000" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // MAX_SIM_BB_EVAL:
  // ----------------
  if ( display_all || NOMAD::string_find ( "MAX_SIM_BB_EVAL ADVANCED \
                                            BLACK-BOXES BLACKBOXES BBEVAL \
                                            MAXIMUM CACHE SIMULATED \
                                            NUMBER EVALUATIONS STOPPING \
                                            TERMINATES TERMINATION" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MAX_SIM_BB_EVAL (advanced)" ) << std::endl
	 << ". maximum number of simulated blackbox evaluations" << std::endl
	 << ". argument: one positive integer" << std::endl
	 << ". no default" << std::endl
	 << ". the same as MAX_BB_EVAL except that it considers" << std::endl
	 << "    initial cache hits (cache points that come from" << std::endl
	 << "    a cache file)" << std::endl
	 << ". simulates the number of blackbox evaluations" << std::endl
	 << "    when no cache file is used" << std::endl
	 << ". example: MAX_SIM_BB_EVAL 1000" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // MAX_TIME:
  // ---------
  if ( display_all || NOMAD::string_find ( "MAX_TIME ADVANCED BLACK-BOXES BLACKBOXES \
                                            MAXIMUM WALL-CLOCK REAL \
                                            STOPPING TERMINATION \
                                            TERMINATES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MAX_TIME (basic)" ) << std::endl
	 << ". maximum wall-clock time in seconds" << std::endl
	 << ". argument: one positive integer" << std::endl
	 << ". no default" << std::endl
	 << ". example: MAX_TIME 3600 # one hour max" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MESH_COARSENING_EXPONENT:
  // -------------------------
  if ( display_all || NOMAD::string_find ( "MESH_COARSENING_EXPONENT ADVANCED \
                                            MADS W+ \\DELTA" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MESH_COARSENING_EXPONENT (advanced)" ) << std::endl
	 << ". mesh coarsening exponent w^+ used to update the mesh" << std::endl
	 << "  after successes (\\Delta^m_{k+1}=\\tau^{w^+}\\Delta^m_k)" << std::endl
	 << ". argument: one nonnegative integer" << std::endl
	 << ". default: 1" << std::endl
	 << ". example: MESH_COARSENING_EXPONENT 0 # the mesh size is" << std::endl
	 << "                                      # not increased" << std::endl
	 << "                                      # after successes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MESH_REFINING_EXPONENT:
  // -----------------------
  if ( display_all || NOMAD::string_find ( "MESH_REFINING_EXPONENT ADVANCED \
                                            MADS W- \\DELTA" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MESH_REFINING_EXPONENT (advanced)" ) << std::endl
	 << ". mesh refining exponent w^- used to update the mesh" << std::endl
	 << "    after failures (\\Delta^m_{k+1} = \\tau^{w^-}\\Delta^m_k)" << std::endl
	 << ". argument: one negative" << std::endl
	 << ". default: -1" << std::endl
	 << ". example: MESH_REFINING_EXPONENT -2" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MESH_UPDATE_BASIS:
  // ------------------
  if ( display_all || NOMAD::string_find ( "MESH_UPDATE_BASIS ADVANCED \
                                            MADS \\TAU \\DELTA" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MESH_UPDATE_BASIS (advanced)" ) << std::endl
	 << ". mesh update basis \\tau used to update the" << std::endl
	 << "    mesh (\\Delta^m_{k+1} = \\tau^w\\Delta^m_k)" << std::endl
	 << ". argument: one positive real" << std::endl
	 << ". default: 4.0" << std::endl
	 << ". example: MESH_UPDATE_BASIS 2.0" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MIN_MESH_SIZE:
  // --------------
  if ( display_all || NOMAD::string_find ( "MIN_MESH_SIZE ADVANCED \
                                            \\DELTA MINIMUM TERMINATION \
                                            STOPPING TERMINATES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MIN_MESH_SIZE (advanced)" ) << std::endl
	 << ". minimum mesh size" << std::endl
	 << ". arguments: same logic as INITIAL_MESH_SIZE" << std::endl
	 << "    (\'r\' can be used)" << std::endl
	 << ". no default" << std::endl
	 << ". example: MIN_MESH_SIZE r1E-5" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MIN_POLL_SIZE:
  // --------------
  if ( display_all || NOMAD::string_find ( "MIN_POLL_SIZE MESH ADVANCED \
                                            \\DELTA^P MINIMUM TERMINATION \
                                            STOPPING TERMINATES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MIN_POLL_SIZE (advanced)" ) << std::endl
	 << ". minimum poll size" << std::endl
	 << ". arguments: same logic as INITIAL_MESH_SIZE" << std::endl
	 << "    (\'r\' can be used)" << std::endl
	 << ". default: 1.0 for integer or binary variables, no default otherwise"
         << std::endl
	 << ". example: MIN_POLL_SIZE r1E-5" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MULTI_F_BOUNDS:
  // ---------------
  if ( display_all || NOMAD::string_find ( "MULTI_F_BOUNDS ADVANCED PARETO \
                                            BI-OBJECTIVES MULTI-OBJECTIVES \
                                            BIOBJECTIVES MULTIOBJECTIVES \
                                            BI-MADS BIMADS SURF" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MULTI_F_BOUNDS (advanced)" ) << std::endl
	 << ". multi-objective optimization: bounds on the two" << std::endl
	 << "    objective functions" << std::endl
	 << ". arguments: 4 reals: f1_min f1_max f2_min f2_max" << std::endl
	 << ". default: none" << std::endl
	 << ". these values are used to display the \'surf\' statistics" << std::endl
	 << "    on Pareto fronts (useful to compare different Pareto" << std::endl
	 << "    fronts)" << std::endl
	 << ". \'surf\' will not be displayed with invalid values (for example"
	 << std::endl
	 << "    if a dominant point has a f2 value greater than f2_max)" << std::endl
	 << ". example: MULTI_F_BOUNDS 0 10 0 10" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MULTI_FORMULATION:
  // -----------------
  if ( display_all || NOMAD::string_find ( "MULTI_FORMULATION ADVANCED PARETO \
                                            BI-OBJECTIVES MULTI-OBJECTIVES\
                                            BIOBJECTIVES MULTIOBJECTIVES \
                                            BI-MADS BIMADS" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MULTI_FORMULATION (advanced)" ) << std::endl
	 << ". multi-objective optimization: single-objective reformulation"
	 << std::endl
	 << ". argument: one string in { \'NORMALIZED\', \'PRODUCT\', \'DIST_L1\',"
	 << std::endl
	 << "                            \'DIST_L2\', \'DIST_LINF\' }" << std::endl
	 << "            (\'NORMALIZED\' and \'DIST_LINF\' are equivalent)" << std::endl
	 << ". default: \'PRODUCT\' or \'DIST_L2\' if VNS_SEARCH is set to \'yes\'"
	 << std::endl
	 << ". example: MULTI_FORMULATION DIST_L1" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }
  
  // MULTI_NB_MADS_RUNS:
  // -------------------
  if ( display_all || NOMAD::string_find ( "MULTI_NB_MADS_RUNS BASIC \
                                            BI-OBJECTIVES MULTI-OBJECTIVES \
                                            BIOBJECTIVES MULTIOBJECTIVES \
                                            BI-MADS BIMADS \
                                            NUMBER" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MULTI_NB_MADS_RUNS (basic)" ) << std::endl
	 << ". multi-objective optimization:" << std::endl
	 << "    number of MADS runs" << std::endl
	 << ". argument: one positive integer" << std::endl
	 << ". default: see user guide" << std::endl
	 << ". example: MULTI_NB_MADS_RUNS 30" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MULTI_OVERALL_BB_EVAL:
  // ---------------------
  if ( display_all || NOMAD::string_find ( "MULTI_OVERALL_BB_EVAL BASIC \
                                            BI-OBJECTIVES MULTI-OBJECTIVES \
                                            BIOBJECTIVES MULTIOBJECTIVES \
                                            BI-MADS BIMADS \
                                            NUMBER BLACK-BOXES BLACKBOXES BBEVAL \
                                            EVALUATIONS TERMINATION \
                                            STOPPING TERMINATES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MULTI_OVERALL_BB_EVAL (basic)" ) << std::endl
	 << ". multi-objective optimization: maximum" << std::endl
	 << "    number of blackbox evaluations" << std::endl
	 << ". argument: one positive integer" << std::endl
	 << ". default: see user guide" << std::endl
	 << ". example: MULTI_OVERALL_BB_EVAL 1000" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // MULTI_USE_DELTA_CRIT:
  // ---------------------
  if ( display_all || NOMAD::string_find ( "MULTI_USE_DELTA_CRITERION ADVANCED PARETO \
                                            BIOBJECTIVES MULTIOBJECTIVES \
                                            BI-MADS BIMADS \
                                            BI-OBJECTIVES MULTI-OBJECTIVES" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "MULTI_USE_DELTA_CRIT (advanced)" ) << std::endl
	 << ". multi-objective optimization: use the delta criterion" << std::endl
	 << "    (can result in a better distributed Pareto front)" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'no\'" << std::endl
	 << ". example: MULTI_USE_DELTA_CRIT yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // OPEN_BRACE:
  // -----------
  if ( display_all || NOMAD::string_find ( "OPEN_BRACES INDENTATION TABULATIONS \
                                            BLOCKS ADVANCED DISPLAY" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "OPEN_BRACE (advanced)" ) << std::endl
	   << ". string displayed at the beginning of indented" << std::endl
	   << "    blocks in full display mode" << std::endl
	   << ". argument: one string" << std::endl
	   << ". default: \'{\'" << std::endl
	   << ". example: OPEN_BRACE Begin" << std::endl
	   << NOMAD::close_block();
      chk = true;
    }

  // OPPORTUNISTIC_CACHE_SEARCH:
  // ---------------------------
  if ( display_all || NOMAD::string_find ( "OPPORTUNISTIC_CACHE_SEARCH ADVANCED \
                                            BLACK-BOXES BLACKBOXES EVALUATIONS \
                                            SUCCESSES CACHE SEARCH" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "OPPORTUNISTIC_CACHE_SEARCH (advanced)" ) << std::endl
	 << ". opportunistic strategy for cache search" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'no\'" << std::endl
	 << ". example: OPPORTUNISTIC_CACHE_SEARCH yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // OPPORTUNISTIC_EVAL:
  // -------------------
  if ( display_all || NOMAD::string_find ( "OPPORTUNISTIC_EVAL BASIC \
                                            BLACK-BOXES BLACKBOXES EVALUATIONS \
                                            SUCCESSES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "OPPORTUNISTIC_EVAL (basic)" ) << std::endl
	 << ". opportunistic strategy (terminate a list of" << std::endl
	 << "    evaluations after successes)" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'yes\'" << std::endl
	 << ". type \'nomad -h opportunistic\' to see advanced options" << std::endl
	 << ". example: OPPORTUNISTIC_EVAL no  # complete evaluations"  << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // OPPORTUNISTIC_LH:
  // -----------------
  if ( display_all || NOMAD::string_find ( "OPPORTUNISTIC_LH BASIC \
                                            BLACK-BOXES BLACKBOXES EVALUATIONS \
                                            SUCCESSES LATIN-HYPERCUBE SEARCH" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "OPPORTUNISTIC_LH (basic)" ) << std::endl
	 << ". opportunistic strategy for Latin-Hypercube search" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: same value as OPPORTUNISTIC_EVAL" << std::endl
	 << ". example: OPPORTUNISTIC_LH no # complete evaluations" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // OPPORTUNISTIC_LUCKY_EVAL:
  // -------------------------
  if ( display_all || NOMAD::string_find ( "OPPORTUNISTIC_LUCKY_EVAL ADVANCED \
                                            BLACK-BOXES BLACKBOXES EVALUATIONS \
                                            SUCCESSES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "OPPORTUNISTIC_LUCKY_EVAL (advanced)" ) << std::endl
	 << ". advanced parameter for the opportunistic" << std::endl
	 << "    strategy (see user guide)" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'no\'" << std::endl
	 << ". example: OPPORTUNISTIC_LUCKY_EVAL yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // OPPORTUNISTIC_MIN_EVAL:
  // -----------------------
  if ( display_all || NOMAD::string_find ( "OPPORTUNISTIC_MIN_EVAL ADVANCED \
                                            BLACK-BOXES BLACKBOXES EVALUATIONS \
                                            SUCCESSES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "OPPORTUNISTIC_MIN_EVAL (advanced)" ) << std::endl
	 << ". advanced parameter for the opportunistic" << std::endl
	 << "    strategy (see user guide)" << std::endl
	 << ". argument: one nonnegative integer" << std::endl
	 << ". no default" << std::endl
	 << ". example: OPPORTUNISTIC_MIN_EVAL 3" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // OPPORTUNISTIC_MIN_F_IMPRVMT:
  // ----------------------------
  if ( display_all || NOMAD::string_find ( "OPPORTUNISTIC_MIN_F_IMPRVMT ADVANCED \
                                            OBJECTIVE \
                                            BLACK-BOXES BLACKBOXES EVALUATIONS \
                                            SUCCESSES IMPROVEMENT" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "OPPORTUNISTIC_MIN_F_IMPRVMT (advanced)" ) << std::endl
	 << ". advanced parameter for the opportunistic" << std::endl
	 << "    strategy (see user guide)" << std::endl
	 << ". argument: one real" << std::endl
	 << ". no default" << std::endl
	 << ". example: OPPORTUNISTIC_MIN_F_IMPRVMT 0.1" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // OPPORTUNISTIC_MIN_NB_SUCCESS:
  // -----------------------------
  if ( display_all || NOMAD::string_find ( "OPPORTUNISTIC_MIN_NB_SUCCESSES ADVANCED \
                                            BLACK-BOXES BLACKBOXES \
                                            EVALUATIONS" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "OPPORTUNISTIC_MIN_NB_SUCCESS (advanced)" ) << std::endl
	 << ". advanced parameter for the opportunistic" << std::endl
	 << "    strategy (see user guide)" << std::endl
	 << ". argument: one nonnegative integer" << std::endl
	 << ". no default" << std::endl
	 << ". example: OPPORTUNISTIC_MIN_NB_SUCCESS 2" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // OPT_ONLY_SGTE:
  // --------------
  if ( display_all || NOMAD::string_find ( "OPT_ONLY_SGTES ADVANCED SURROGATES \
                                            BLACK-BOXES	BLACKBOXES \
                                            SGTES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "OPT_ONLY_SGTE (advanced)" ) << std::endl
	 << ". NOMAD will only minimize the surrogate" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". SGTE_EXE or HAS_SGTE must be defined" << std::endl
	 << ". default: \'no\'" << std::endl
	 << ". example: OPT_ONLY_SGTE yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // PERIODIC_VARIABLE:
  // ------------------
  if ( display_all || NOMAD::string_find ( "PERIODIC_VARIABLE ADVANCED \
                                            BOUNDS LB UB CYCLIC MADS" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "PERIODIC_VARIABLE (advanced)" ) << std::endl
	   << ". specify that some variables are periodic" << std::endl
	   << ". arguments: variable indexes" << std::endl
	   << ". no default" << std::endl
	   << ". bounds must be defined for these variables" << std::endl
	   << ". " << NOMAD::open_block ( "examples" )
	   << "PERIODIC_VARIABLE *   # all variables are periodic" << std::endl
	   << "PERIODIC_VARIABLE 0-1 # 2 first var. are periodic" << std::endl
	   << NOMAD::close_block()
	   << NOMAD::close_block();
      chk = true;
    }

  // POINT_DISPLAY_LIMIT:
  // --------------------
  if ( display_all || NOMAD::string_find ( "POINT_DISPLAY_LIMIT OUTPUTS \
                                            ADVANCED PRECISION" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "POINT_DISPLAY_LIMIT (advanced)" ) << std::endl
	 << ". maximum number of point coordinates" << std::endl
	 << "    that are displayed" << std::endl
	 << ". argument: one positive integer" << std::endl
	 << ". default: 20" << std::endl
	 << ". example: POINT_DISPLAY_LIMIT 10" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // RHO:
  // ----
  if ( display_all || NOMAD::string_find ( "RHO ADVANCED MADS CONSTRAINTS \
                                            PROGRESSIVE-BARRIER PB PEB	\
                                            FILTER TRIGGER" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "RHO (advanced)" ) << std::endl
	 << ". rho parameter of the progressive barrier" << std::endl
	 << ". argument: one nonnegative real" << std::endl
	 << ". default: 0.1" << std::endl
	 << ". example: RHO 0.5" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // SCALING:
  // --------
  if ( display_all || NOMAD::string_find ( "SCALING SCALE ADVANCED \
                                            FILES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "SCALING (advanced)" ) << std::endl
	 << ". variable scaling" << std::endl
	 << ". arguments: variable indexes and values" << std::endl
	 << ". no default" << std::endl
	 << ". variables are multiplied by these values: they are scaled" << std::endl
	 << "    before an evaluation and the call to Evaluator::eval_x()," << std::endl
	 << "    and unscaled after the evaluation" << std::endl
	 << ". all NOMAD outputs (including files) display unscaled values" << std::endl
	 << ". all variable-related parameters (bounds, starting points," << std::endl
	 << "    fixed variables) must be specified without scaling" << std::endl
	 << ". can be given by a text file containing DIMENSION entries" << std::endl
	 << "    (use \'-\' for unscaled variables)" << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "SCALING ( 0.1 - 100 ) # variables 0 and 2 are scaled" << std::endl
	 << "SCALING scaling.txt   # with a file" << std::endl
	 << "SCALING 0-1 10.0      # 2 first variables scaled" << std::endl
	 << "                      # by factor 10" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }

  // SEC_POLL_DIR_TYPES:
  // -------------------
  if ( display_all || NOMAD::string_find ( "SEC_POLL_DIR_TYPES ADVANCED MADS \
                                            POLL DIRECTIONS PB PEB \
                                            PROGRESSIVE-BARRIER" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "SEC_POLL_DIR_TYPES (advanced)" ) << std::endl
	 << ". types of directions for the secondary poll" << std::endl
	 << ". arguments: same logic as DIRECTION_TYPE" << std::endl
	 << ". default: see user guide" << std::endl
	 << ". example: SEC_POLL_DIR_TYPES ORTHO 1" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // SEED:
  // -----
  if ( display_all || NOMAD::string_find ( "SEED BASIC \
                                            RANDOM FILES LT-MADS LTMADS \
                                            LATIN-HYPERCUBE \
                                            SAMPLING SEARCH" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "SEED (basic)" ) << std::endl
	 << ". random seed" << std::endl
	 << ". argument: one integer or the string \'NONE\'" << std::endl
	 << ". default: \'NONE\'" << std::endl
	 << ". the random seed is different at each run if" << std::endl
	 << "    \'NONE\' or a negative integer is entered" << std::endl
	 << ". the seed is used in the output file names" << std::endl
	 << ". randomness is used for LT-MADS directions" << std::endl
	 << "    and for Latin-Hypercube search" << std::endl
	 << ". example: SEED 100" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // SGTE_CACHE_FILE:
  // ----------------
  if ( display_all || NOMAD::string_find ( "SGTE_CACHE_FILE ADVANCED \
                                            SURROGATES SGTES \
                                            FILES OUTPUTS" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "SGTE_CACHE_FILE (advanced)" ) << std::endl
	 << ". surrogate cache file; cannot be the same" << std::endl
	 << "    as CACHE_FILE" << std::endl
	 << ". argument: one string" << std::endl
	 << ". no default" << std::endl
	 << ". points already in the file will be tagged" << std::endl
	 << "    as surrogate evaluations" << std::endl
	 << ". example: SGTE_CACHE_FILE sgte_cache.bin" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // SGTE_COST:
  // ----------
  if ( display_all || NOMAD::string_find ( "SGTE_COST SURROGATES SGTES ADVANCED \
                                            BLACK-BOXES BLACKBOXES" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "SGTE_COST (advanced)" ) << std::endl
	   << ". cost of the surrogate function relatively" << std::endl
	   << "    to the true function" << std::endl
	   << ". argument: one nonnegative integer" << std::endl
	   << ". default: infinity (no surrogate cost)" << std::endl
	   << ". " << NOMAD::open_block ( "examples" )
	   << "SGTE_COST 3   # three surrogate evaluations" << std::endl
	   << "              # count as one blackbox" << std::endl
	   << "              # evaluation (the surrogate" << std::endl
	   << "              # is three times faster)" << std::endl
	   << "SGTE_COST -1  # set to infinity" << std::endl
	   << NOMAD::close_block()
	   << NOMAD::close_block();
      chk = true;
    }

  // SGTE_EVAL_SORT:
  // ---------------
  if ( display_all || NOMAD::string_find ( "SGTE_EVAL_SORT ADVANCED SURROGATES \
                                            SGTES BLACK-BOXES \
                                            BLACKBOXES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "SGTE_EVAL_SORT (advanced)" ) << std::endl
	 << ". if surrogate is used to sort the trial points" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'yes\'" << std::endl
	 << ". example: SGTE_EVAL_SORT NO" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // SGTE_EXE:
  // ----------
  if ( display_all || NOMAD::string_find ( "SGTE_EXE HAS_SGTE ADVANCED SURROGATES \
                                            SGTES BLACK-BOXES	BLACKBOXES \
                                            FILES EXECUTABLE" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "SGTE_EXE (advanced)" ) << std::endl
	 << ". to indicate a surrogate executable" << std::endl
	 << ". arguments: one or two strings" << std::endl
	 << ". no default" << std::endl
	 << ". surrogate(s) and blackbox(es) must have the same" << std::endl
	 << "    number of outputs" << std::endl
	 << ". if surrogates are used, every blackbox executable" << std::endl
	 << "    must have a surrogate"
	 << std::endl
	 << ". automatically sets HAS_SGTE to \'yes\'" << std::endl
	 << ". " << NOMAD::open_block ( "examples" ) << std::endl
	 << "SGTE_EXE b1.exe s1.exe # \'s1.exe\' is a surrogate" << std::endl
	 << "                       # for \'b1.exe\'" << std::endl << std::endl
	 << "SGTE_EXE sgte.exe      # only if one blackbox" << std::endl
	 << "                       # executable is used" << std::endl 
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }

  // SNAP_TO_BOUNDS:
  // ---------------
  if ( display_all || NOMAD::string_find ( "SNAP_TO_BOUNDS ADVANCED" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "SNAP_TO_BOUNDS (advanced)" ) << std::endl
	   << ". if \'yes\', snap to bounds points generated" << std::endl
	   << "    outside boundaries" << std::endl
	   << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	   << ". default: \'yes\'" << std::endl
	   << ". example: SNAP_TO_BOUNDS no" << std::endl
	   << NOMAD::close_block();
      chk = true;
    }

  // SOLUTION_FILE:
  // --------------
  if ( display_all || NOMAD::string_find ( "SOLUTION_FILE BASIC \
                                            FILES OUTPUTS" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "SOLUTION_FILE (basic)" ) << std::endl
	 << ". file containing the solution" << std::endl
	 << ". argument: one string" << std::endl
	 << ". no default" << std::endl
	 << ". the seed is added to the file name if" << std::endl
	 << "  ADD_SEED_TO_FILE_NAMES=\'yes\'" << std::endl
	 << ". example: SOLUTION_FILE sol.txt" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // SPECULATIVE_SEARCH:
  // -------------------
  if ( display_all || NOMAD::string_find ( "SPECULATIVE_SEARCH MADS OPTIMISTIC \
                                            ADVANCED" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "SPECULATIVE_SEARCH (advanced)" ) << std::endl
	 << ". MADS speculative_search (optimistic strategy)" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'yes\'" << std::endl
	 << ". example: SPECULATIVE_SEARCH no" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // STAT_SUM_TARGET:
  // ----------------
  if ( display_all || NOMAD::string_find ( "STAT_SUM_TARGET ADVANCED TERMINATION \
                                            STOPPING TERMINATES STATS" , param_names ) )
    {
      _out << std::endl
	   << NOMAD::open_block ( "STAT_SUM_TARGET (advanced)" ) << std::endl
	   << ". MADS terminates if STAT_SUM reaches the value of this" << std::endl
	   << "    parameter (STAT_SUM is one of the possible outputs" << std::endl
	   << "    defined in BB_OUTPUT_TYPE)" << std::endl
	   << ". argument: one real" << std::endl
	   << ". no default" << std::endl
	   << ". example: STAT_SUM_TARGET 100.0" << std::endl
	   << NOMAD::close_block();
      chk = true;
    }
  
  // STATS_FILE:
  // -----------
  if ( display_all || NOMAD::string_find ( "STATS_FILE BASIC \
                                            FILES OUTPUTS DISPLAY_STATS" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "STATS_FILE (basic)" ) << std::endl
	 << ". file containing all successes with the same" << std::endl
	 << "    format than DISPLAY_STATS" << std::endl
	 << ". arguments: one string (file name) and one" << std::endl
	 << "    list of strings (stats)" << std::endl
	 << ". no default" << std::endl
	 << ". the seed is added to the file name if" << std::endl
	 << "    ADD_SEED_TO_FILE_NAMES=\'yes\'" << std::endl
	 << ". example: STATS_FILE log.txt BBE SOL f=%.2EOBJ" << std::endl
	<< NOMAD::close_block();
    chk = true;
  }

  // STOP_IF_FEASIBLE:
  // -----------------
  if ( display_all || NOMAD::string_find ( "STOP_IF_FEASIBLE ADVANCED \
                                            TERMINATES TERMINATION STOPPING" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "STOP_IF_FEASIBLE (advanced)" ) << std::endl
	 << ". the algorithm terminates if it generates" << std::endl
	 << "    a feasible solution" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'no\'" << std::endl
	 << ". example: STOP_IF_FEASIBLE yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // TMP_DIR:
  // --------
  if ( display_all || NOMAD::string_find ( "TMP_DIR BASIC PATH TEMPORARY \
                                            DIRECTORY FILES \
                                            BLACK-BOXES BLACKBOXES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "TMP_DIR (basic)" ) << std::endl
	 << ". temporary directory for blackbox input/output files"  << std::endl
	 << ". argument: one string indicating a directory" << std::endl
	 << ". default: problem directory" << std::endl
	 << ". improved performance with a local temporary directory" << std::endl
	 << ". example: TMP_DIR /tmp" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // UNDEF_STR:
  // ----------
  if ( display_all || NOMAD::string_find ( "UNDEF_STR ADVANCED \
                                            UNDEFINED DISPLAY REALS" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "UNDEF_STR (advanced)" ) << std::endl
	 << ". string used to display undefined real values" << std::endl
	 << ". argument: one string" << std::endl
	 << ". default: \"-\"" << std::endl
	 << ". example: UNDEF_STR Undefined" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // UPPER_BOUND:
  // ------------
  if ( display_all || NOMAD::string_find ( "UPPER_BOUND BASIC UB BOUNDS" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "UPPER_BOUND (basic)" ) << std::endl
	 << ". upper bounds for each variable" << std::endl
	 << ". same logic as parameter LOWER_BOUND" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // USER_CALLS_ENABLED:
  // -------------------
  if ( display_all || NOMAD::string_find ( "USER_CALLS_ENABLED USER_CALLS_DISABLED \
                                            ADVANCED LIBRARY" ,
					   param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "USER_CALLS_ENABLED (advanced)" ) << std::endl
	 << ". if \'no\', the automatic calls to user" << std::endl
	 << "    functions are disabled" << std::endl
	 << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	 << ". default: \'yes\'" << std::endl
	 << ". example: USER_CALLS_ENABLED yes" << std::endl
	 << NOMAD::close_block();
    chk = true;
  }

  // VARIABLE_GROUP:
  // --------------
  if ( display_all || NOMAD::string_find ( "VARIABLE_GROUP GROUPS PSD-MADS PSDMADS \
                                            ADVANCED" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "VARIABLE_GROUP (advanced)" ) << std::endl
	 << ". defines groups of variables" << std::endl
	 << ". MADS directions are applied separately for" << std::endl
	 << "    each group" << std::endl
	 << ". also used by PSD-MADS (not yet implemented)" << std::endl
	 << ". groups cannot include fixed variables" << std::endl
	 << ". arguments: variable indexes" << std::endl
	 << ". default groups are created for different types" << std::endl
	 << "    of variables" << std::endl
	 << ". no other default" << std::endl
	 << ". advanced options only available in library mode" << std::endl
	 << "    (see user guide)" << std::endl
	 << ". " << NOMAD::open_block ( "examples" )
	 << "VARIABLE_GROUP 2-5" << std::endl
	 << "VARIABLE_GROUP 0 1 3" << std::endl
      	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }

  // VNS_SEARCH:
  // -----------
  if ( display_all || NOMAD::string_find ( "VNS_SEARCH VARIABLE NEIGHBORHOOD \
                                            METAHEURISTICS META-HEURISTICS \
                                            GLOBAL BASIC \
                                            TRIGGER" ,
					   param_names ) ) {

    if ( !NOMAD::string_find ( "RHO" , param_names ) ) {

      _out << std::endl
	   << NOMAD::open_block ( "VNS_SEARCH (basic)" ) << std::endl
	   << ". Variable Neighborhood Search (VNS) search" << std::endl
	   << ". argument: one boolean (\'yes\' or \'no\')" << std::endl
	   << "         or one real in [0;1] for the VNS trigger" << std::endl
	   << ". default: \'no\' (same as 0.0)" << std::endl
	   << ". the VNS trigger is the maximum desired ratio of" << std::endl
	   << "    VNS blackbox evaluations over the total number" << std::endl
	   << "    of blackbox evaluations" << std::endl
	   << ". the VNS search is never executed with a null trigger" << std::endl
	   << "    while a value of 1 allows the search at every" << std::endl
	   << "    iteration" << std::endl
	   << ". if VNS_SEARCH=\'yes\', the default value of 0.75 is" << std::endl
	   << "    taken for the trigger" << std::endl
	   << ". VNS search uses the surrogate if HAS_SGTE or" << std::endl
	   << "    SGTE_EXE is defined" << std::endl
	   << ". " << NOMAD::open_block ( "examples" )
	   << "VNS_SEARCH yes  # VNS trigger of 75%" << std::endl
	   << "VNS_SEARCH 0.5  # VNS trigger of 50%" << std::endl
	   << NOMAD::close_block()
	   << NOMAD::close_block();
      chk = true;
    }
  }
  
  // X0:
  // ---
  if ( display_all || NOMAD::string_find ( "X0 STARTING POINT BASIC \
                                            VARIABLES \
                                            LH LATIN-HYPERCUBE \
                                            CACHE FILES" , param_names ) ) {
    _out << std::endl
	 << NOMAD::open_block ( "X0 (basic)" ) << std::endl
	 << ". starting point(s)" << std::endl
	 << ". arguments: text file name," << std::endl
	 << "          or cache file name," << std::endl
	 << "          or DIMENSION reals," << std::endl
	 << "          or indexed values" << std::endl
	 << ". default: best point from a cache file or from" << std::endl
	 << "           an initial LH search" << std::endl
	 << ". do not use a surrogate cache file" << std::endl
	 << "    (even if OPT_ONLY_SGTE=\'yes\')" << std::endl
	 << ". more than one starting point can be defined (all points" << std::endl
	 << "    are evaluated: x0 evaluations are not opportunistic)" << std:: endl
	 << ". a text file can describe more than one point" << std::endl
	 << ". may be infeasible, but can only violate PB, F, or PEB" << std::endl
	 << "    constraints" << std::endl
	 << ". cannot be outside bounds" << std::endl
	 << ". must respect fixed variables (param. FIXED_VARIABLE)" << std::endl
	 << ". " << NOMAD::open_block ("examples") << std::endl
	 << "X0 x0.txt     # text file with a multiple" << std::endl
	 << "              # of DIMENSION values" << std::endl << std::endl
	 << "X0   * 0.0    # first starting point" << std::endl
	 << "X0 1 * 1.0    # second starting point" << std::endl << std::endl
	 << "X0 ( 0 1 2 )  # if DIMENSION=3" << std::endl << std::endl
	 << "see other examples in user guide" << std::endl
	 << NOMAD::close_block()
	 << NOMAD::close_block();
    chk = true;
  }
  
  // last display:
  if ( !chk ) {

    std::string pname = ( pnames.size() == 1 ) ?
      ("\'" + *pnames.begin() + "\'") :
      "the specified list of parameter names";

    _out << std::endl << "no help available for " << pname << std::endl
	 << "help example: \'nomad -h mesh\' displays help on the mesh parameters"
	 << std::endl;
  }
 
  _out.close_block();
}
