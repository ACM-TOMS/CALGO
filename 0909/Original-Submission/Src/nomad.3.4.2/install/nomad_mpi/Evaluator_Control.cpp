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
  \file   Evaluator_Control.cpp
  \brief  Control of the blackbox evaluations (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-15
  \see    Evaluator_Control.hpp
*/
#include "Evaluator_Control.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
bool NOMAD::Evaluator_Control::_force_quit = false;

/*---------------------------------------------------------*/
/*                       constructor                       */
/*---------------------------------------------------------*/
NOMAD::Evaluator_Control::Evaluator_Control
( const NOMAD::Parameters & p          ,
  NOMAD::Stats            & stats      ,
  NOMAD::Evaluator        * ev         ,   // can be NULL
  NOMAD::Cache            * cache      ,   // can be NULL
  NOMAD::Cache            * sgte_cache   ) // can be NULL
  : _p                ( p          ) ,
    _ev               ( ev         ) ,
    _del_ev           ( false      ) ,
    _cache            ( cache      ) ,
    _sgte_cache       ( sgte_cache ) ,
    _del_cache        ( false      ) ,
    _del_sgte_cache   ( false      ) ,
#ifdef USE_MPI
    _eval_in_progress ( NULL       ) ,
    _nb_in_progress   ( 0          ) ,
    _elop_tag         ( 0          ) ,
    _slaves_elop_tags ( NULL       ) ,
    _slave            ( NULL       ) ,
#endif
    _stats            ( stats      ) ,
    _last_stats_tag   ( -1         ) ,
    _last_stats_bbe   ( -1         )
{
  NOMAD::Evaluator_Control::_force_quit = false;

  // Evaluator init:
  if ( !_ev ) {
    _ev = ( _p.get_index_obj().size() > 1 ) ? new NOMAD::Multi_Obj_Evaluator ( p ):
                                              new NOMAD::Evaluator           ( p );
    _del_ev = true;
  }

  if ( NOMAD::Slave::is_master() ) {

#ifdef USE_MPI

    int np = NOMAD::Slave::get_nb_processes();

    _eval_in_progress = new NOMAD::Eval_Point * [np];
    _slaves_elop_tags = new int                 [np];
    for ( int i = 0 ; i < np ; ++i ) {
      _eval_in_progress[i] = NULL;
      _slaves_elop_tags[i] = -1;
    }
   
    _slave = new NOMAD::Slave ( _p , _ev );

#endif
    
    const NOMAD::Display & out = _p.out();

    // caches creation:
    if ( !_cache ) {
      _cache     = new NOMAD::Cache ( out , NOMAD::TRUTH );
      _del_cache = true;
    }
    if ( !_sgte_cache ) {
      _sgte_cache     = new NOMAD::Cache ( out , NOMAD::SGTE );
      _del_sgte_cache = true;
    }

    // caches init (we only load cache file points with m blackbox outputs):
    std::string    file_name;
    int            m              = p.get_bb_nb_outputs();
    NOMAD::dd_type display_degree = out.get_gen_dd();
    
    if ( !_p.get_cache_file().empty() ) {
      file_name = _p.get_problem_dir() + _p.get_cache_file();
      if ( !_cache->load ( file_name , &m , display_degree == NOMAD::FULL_DISPLAY )
	   && display_degree != NOMAD::NO_DISPLAY )
	out << std::endl
	    << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
	    << "): could not load (or create) the cache file " << file_name
	    << std::endl << std::endl;
    }
    
    if ( !_p.get_sgte_cache_file().empty() ) {
      file_name = _p.get_problem_dir() + _p.get_sgte_cache_file();
      if ( !_sgte_cache->load ( file_name ,
				&m        ,
				display_degree==NOMAD::FULL_DISPLAY ) &&
	   display_degree != NOMAD::NO_DISPLAY )
	out << std::endl << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
	    << "): could not load (or create) the surrogate cache file "
	    << file_name << std::endl << std::endl;
    }
  }
}

/*---------------------------------------------------------*/
/*                        destructor                       */
/*---------------------------------------------------------*/
NOMAD::Evaluator_Control::~Evaluator_Control ( void )
{
  if ( _del_ev )
    delete _ev;
  if ( _del_cache )
    delete _cache;
  if ( _del_sgte_cache )
    delete _sgte_cache;

  clear_eval_lop();

#ifdef USE_MPI

  if ( _eval_in_progress ) {
    int np = NOMAD::Slave::get_nb_processes();
    for ( int i = 0 ; i < np ; ++i )
      if ( _eval_in_progress[i] && !_eval_in_progress[i]->is_in_cache() )
	delete _eval_in_progress[i];
    delete [] _eval_in_progress;
  }
  if ( _slaves_elop_tags )
    delete [] _slaves_elop_tags;

  delete _slave;

#endif
}

/*---------------------------------------------------------*/
/*                     save the caches                     */
/*---------------------------------------------------------*/
bool NOMAD::Evaluator_Control::save_caches ( bool overwrite )
{
  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_gen_dd();

  bool b1 = _cache->save      ( overwrite , display_degree == NOMAD::FULL_DISPLAY );
  bool b2 = _sgte_cache->save ( overwrite , display_degree == NOMAD::FULL_DISPLAY );

  if ( !b1 && display_degree != NOMAD::NO_DISPLAY )
    out << std::endl << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
	<< "): could not save the cache file "
	<< _p.get_problem_dir() << _p.get_cache_file()
	<< std::endl << std::endl;

  if ( !b2 && display_degree != NOMAD::NO_DISPLAY )
    out << std::endl
	<< "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
	<< "): could not save the surrogate cache file "
	<< _p.get_problem_dir() << _p.get_sgte_cache_file()
	<< std::endl << std::endl; 
  return b1 && b2;
}

/*---------------------------------------------------------*/
/*    process an already evaluated Eval_Point (private)    */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::process_eval_point
( const NOMAD::Eval_Point & x            ,
  NOMAD::Barrier          & barrier      ,
  NOMAD::Pareto_Front     * pareto_front ) const
{
  // insertion of the Eval_Point in the barriers:
  barrier.insert ( x );

  if ( x.get_eval_type() == NOMAD::TRUTH || _p.get_opt_only_sgte() ) {

    // multi-objective:
    if ( pareto_front ) {

      // insertion of the Eval_Point in the Pareto front:
      if ( x.is_feasible ( _p.get_h_min() ) &&
	   pareto_front->insert ( x )       &&
	   _p.get_user_calls_enabled()         )
	_ev->update_success ( _stats , x );

    }
   
    // single-objective: call virtual method Evaluator::update_success():
    else if ( _p.get_user_calls_enabled() &&
	      barrier.get_one_eval_succ() == NOMAD::FULL_SUCCESS )
      _ev->update_success ( _stats , x );
  } 
}

/*---------------------------------------------------------*/
/*  update barrier b1 from points in barrier b2 and treat  */
/*  these points as evaluations (used in VNS search)       */
/*---------------------------------------------------------*/
NOMAD::success_type NOMAD::Evaluator_Control::process_barrier_points
( NOMAD::Barrier       & b1             ,
  const NOMAD::Barrier & b2             ,
  NOMAD::Pareto_Front  * pareto_front   ,
  NOMAD::dd_type         display_degree ,
  NOMAD::search_type     search           ) const
{
  b1.reset_success();

  NOMAD::Eval_Point                       * modifiable_x;
  NOMAD::success_type                       one_eval_succ;
  bool                                      opt_only_sgte = _p.get_opt_only_sgte();
  const std::string                       & his_file      = _p.get_history_file();
  const NOMAD::Eval_Point                 * last_success  = NULL;
  const std::list<const NOMAD::Eval_Point *>
                                          & all_inserted  = b2.get_all_inserted();
  std::list<const NOMAD::Eval_Point *>::const_iterator
                                            it , end = all_inserted.end();
  for ( it = all_inserted.begin() ; it != end ; ++it ) {

    // insertion in barrier:
    modifiable_x = &NOMAD::Cache::get_modifiable_point ( **it );

    modifiable_x->set_direction          ( NULL                              );
    modifiable_x->set_mesh_index         ( NULL                              );
    modifiable_x->set_poll_center_type   ( NOMAD::UNDEFINED_POLL_CENTER_TYPE );
    modifiable_x->set_user_eval_priority ( NOMAD::Double()                   );
    modifiable_x->set_rand_eval_priority ( NOMAD::Double()                   );

    // process evaluation point:
    process_eval_point ( **it , b1 , pareto_front );
    
    one_eval_succ = b1.get_one_eval_succ();
    if ( one_eval_succ != NOMAD::UNSUCCESSFUL && one_eval_succ >= b1.get_success() )
      last_success = *it;
    
    // update the history file:
    if ( !his_file.empty() && ( opt_only_sgte ||
				(*it)->get_eval_type() == NOMAD::TRUTH ) )
      write_sol_or_his_file ( _p.get_problem_dir() + his_file , **it , false );
  }

  NOMAD::success_type success = b1.get_success();

  // display and save only the last success:
  if ( last_success )
    display_eval_result ( *last_success  ,
			  display_degree ,
			  search         ,
			  success        ,
			  success          );
  
  // barrier update:
  b1.update_and_reset_success();

  return success;
}

/*---------------------------------------------------------*/
/*      count the output stats (STAT_SUM and STAT_AVG)     */
/*      (private)                                          */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::count_output_stats ( const NOMAD::Eval_Point & x )
{
  const NOMAD::Point & bbo   = x.get_bb_outputs();
  int                  i_sum = _p.get_index_stat_sum();
  int                  i_avg = _p.get_index_stat_avg();

  // STAT_SUM:
  if ( i_sum >= 0 )
    _stats.update_stat_sum ( bbo[i_sum] );

  // STAT_AVG:
  if ( i_avg >= 0 )
    _stats.update_stat_avg ( bbo[i_avg] );
}

/*-------------------------------------------------------------------*/
/*                file displays for parameter STATS_FILE             */
/*-------------------------------------------------------------------*/
void NOMAD::Evaluator_Control::stats_file ( const std::string       & file_name ,
					    const NOMAD::Eval_Point * x         ,
					    const NOMAD::Point      * multi_obj ) const
{
  std::string   fn = _p.get_problem_dir() + file_name;
  std::ofstream fout ( fn.c_str() , std::ios::app );

  if ( !fout.fail() ) {
    fout.setf      ( std::ios::fixed             );
    fout.precision ( NOMAD::DISPLAY_PRECISION_BB );
    display_stats ( false , fout , _p.get_stats_file() , x , multi_obj );
  }
  else {
    const NOMAD::Display & out = _p.out();
    if ( out.get_gen_dd() != NOMAD::NO_DISPLAY )
      out << std::endl
	  << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
	  << "): could not save information in stats file \'"
	  << file_name << "\'" << std::endl << std::endl;
  }
  fout.close();
}

/*-------------------------------------------------------------------*/
/*  display stats during Mads::run() for minimal and normal display  */
/*-------------------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_stats
( bool                           header    ,
  const NOMAD::Display         & out       ,
  const std::list<std::string> & stats     ,
  const NOMAD::Eval_Point      * x         ,
  const NOMAD::Point           * multi_obj   ) const
{
  if ( stats.empty() ) {
    out << std::endl;
    return;
  }

  if ( header )
    out << std::endl;

  NOMAD::Double            f;
  const NOMAD::Point     * sol       = NULL;
  const NOMAD::Point     * bbo       = NULL;
  const NOMAD::Signature * signature = NULL;
  int                      bbe       = _stats.get_bb_eval();
  int                      i;

  // this integer is used for the default width display
  // of the various stats on the number of evaluations:
  int max_bbe = _p.get_max_bb_eval();
  if ( _p.get_max_sgte_eval() > max_bbe )
    max_bbe = _p.get_max_sgte_eval();
  if ( _p.get_max_sim_bb_eval() > max_bbe )
    max_bbe = _p.get_max_sim_bb_eval();
  if ( _p.get_max_eval() > max_bbe )
    max_bbe = _p.get_max_eval();

  if ( x ) {
    signature       = x->get_signature();
    f               = x->get_f();
    sol             = x;
    bbo             = &(x->get_bb_outputs());
    _last_stats_tag = x->get_tag();
    _last_stats_bbe = bbe;
  }

  std::string s1 , format;
  std::list<std::string>::const_iterator it , end = stats.end();
  for ( it = stats.begin() ; it != end ; ++it ) {

    if ( it->empty() )
      out << "\t";
    
    else {

      if ( header ) {
	s1 = *it;
	NOMAD::Display::extract_display_format ( s1 , format );
	out << s1;
      }

      else {
	
	// get the stats type:
	NOMAD::display_stats_type dst
	  = NOMAD::Display::get_display_stats_type ( *it );
	
	// some stats types are disables in the multi-objective case:
	if ( multi_obj &&
	     ( dst == NOMAD::DS_SIM_BBE  ||
	       dst == NOMAD::DS_BBE      ||
	       dst == NOMAD::DS_SGTE     ||
	       dst == NOMAD::DS_EVAL     ||
	       dst == NOMAD::DS_TIME     ||
	       dst == NOMAD::DS_STAT_SUM ||
	       dst == NOMAD::DS_STAT_AVG    ) )
	  dst = NOMAD::DS_UNDEFINED;

	// display the stats:
	switch ( dst ) {
	case NOMAD::DS_UNDEFINED:
	  s1 = *it;
	  NOMAD::Display::extract_display_format ( s1 , format );
	  out << s1;
	  break;
	case NOMAD::DS_OBJ:
	  if ( multi_obj )
	    display_stats_point ( out , stats , it , multi_obj );
	  else {
	    display_stats_real ( out , f , format );
	    format.clear();
	  }
	  break;
	case NOMAD::DS_MESH_INDEX:
	  display_stats_int ( out                           ,
			      NOMAD::Mesh::get_mesh_index() ,
			      10*L_LIMITS                   ,
			      format                          );
	  format.clear();
	  break;
	case NOMAD::DS_DELTA_M:
	case NOMAD::DS_MESH_SIZE:
	  {
	    if ( signature ) {
	      NOMAD::Point delta_m;
	      signature->get_mesh().get_delta_m ( delta_m ,
						  NOMAD::Mesh::get_mesh_index() );
	      display_stats_point ( out , stats , it , &delta_m );
	    }
	    else
	      out << "-";
	  }
	  break;
	case NOMAD::DS_DELTA_P:
	case NOMAD::DS_POLL_SIZE:
	  {
	    if ( signature ) {
	      NOMAD::Point delta_p;
	      signature->get_mesh().get_delta_p ( delta_p ,
						  NOMAD::Mesh::get_mesh_index() );
	      display_stats_point ( out , stats , it , &delta_p );
	    }
	    else
	      out << "-";
	  }
	  break;
	case NOMAD::DS_SIM_BBE:
	  display_stats_int ( out , _stats.get_sim_bb_eval() , max_bbe , format );
	  format.clear();
	  break;
	case NOMAD::DS_BBE:
	  display_stats_int ( out , bbe , max_bbe , format );
	  format.clear();
	  break;
	case NOMAD::DS_SGTE:
	  display_stats_int ( out , _stats.get_sgte_eval() , max_bbe , format );
	  format.clear();
	  break;
	case NOMAD::DS_EVAL:
	  display_stats_int ( out , _stats.get_eval() , max_bbe , format );
	  format.clear();
	  break;
	case NOMAD::DS_TIME:
	  display_stats_int ( out , _stats.get_real_time() , 3600 , format );
	  format.clear();
	  break;
	case NOMAD::DS_STAT_SUM:
	  display_stats_real ( out , _stats.get_stat_sum() , format );
	  format.clear();
	  break;
	case NOMAD::DS_STAT_AVG:
	  display_stats_real ( out , _stats.get_stat_avg() , format );
	  format.clear();
	  break;
	case NOMAD::DS_BBO:
	  display_stats_point ( out , stats , it , bbo );
	  break;
	case NOMAD::DS_SOL:
	  display_stats_point ( out , stats , it , sol );
	  break;
	case NOMAD::DS_VAR:
	  ++it;
	  NOMAD::atoi ( *it , i );
	  if ( sol )
	    display_stats_real ( out , (*sol)[i] , format );
	  else
	    out << "-";
	  format.clear();
	  break;
	}
      }
    }
  }
  
  if ( !header )
    out << std::endl;
}

/*-----------------------------------------------------*/
/*  display a real with DISPLAY_STATS (or STATS_FILE)  */
/*-----------------------------------------------------*/
void NOMAD::Evaluator_Control::display_stats_real
( const NOMAD::Display & out    ,
  const NOMAD::Double  & d      ,
  const std::string    & format   ) const
{
  if ( format.empty() ) {
    std::string format2 = "%0." + NOMAD::itos(DISPLAY_PRECISION_STD) + "g";
    d.display ( out , format2 );
  }
  else
    d.display ( out , format );
}

/*---------------------------------------------------------*/
/*  display an integer with DISPLAY_STATS (or STATS_FILE)  */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_stats_int
( const NOMAD::Display & out    ,
  int                    i      ,
  int                    max_i  ,
  const std::string    & format   ) const
{
  if ( format.empty() )
    out.display_int_w ( i , max_i );
  else {
    NOMAD::Double d = i;
    d.display ( out , format );
  }
}

/*---------------------------------------------------------*/
/*    display a point with DISPLAY_STATS (or STATS_FILE)   */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_stats_point 
( const NOMAD::Display                   & out           ,
  const std::list<std::string>           & display_stats ,
  std::list<std::string>::const_iterator & it            ,
  const NOMAD::Point                     * x               ) const
{
  if ( x ) {

    int n = x->size();

    // s1 is the string displayed befores and after
    // one coordinate (it may include format):
    std::string s1;
    if ( it != display_stats.begin() ) {
      s1 = *(--it);
      ++it;
    }

    // extract the display format from s1:
    std::string format;
    if ( !s1.empty() )
      NOMAD::Display::extract_display_format ( s1 , format );

    // s2 is the string displayed between two coordinates:
    std::string s2;
    ++it;
    if ( it != display_stats.end() )
      s2 = *it;
    else if ( s2.empty() )
      --it;
    for ( int i = 0 ; i < n ; ++i ) {
      if ( !s1.empty() && i > 0 )
	out << s1;

      display_stats_real ( out , (*x)[i] , format );

      if ( !s1.empty() )
	out << s1;
      if ( !s2.empty() && i < n-1 )
	out << " " << s2;
      out << " ";
    }
  }
}

/*------------------------------------------*/
/*  save the solution file (SOLUTION_FILE)  */
/*------------------------------------------*/
void NOMAD::Evaluator_Control::write_solution_file
( const NOMAD::Eval_Point & x ) const
{
  const std::string & sol_file = _p.get_solution_file();
  if ( !sol_file.empty() && x.is_feasible ( _p.get_h_min() ) )
    write_sol_or_his_file ( _p.get_problem_dir() + sol_file , x , true );
}

/*----------------------------------------------*/
/*     save the solution file  (SOLUTION_FILE)  */
/*  or update the history file (HISTORY_FILE )  */
/*  (private)                                   */
/*----------------------------------------------*/
void NOMAD::Evaluator_Control::write_sol_or_his_file 
( const std::string       & file_name ,
  const NOMAD::Eval_Point & x         ,
  bool                      is_sol      ) const
{
  // if is_sol == true: save the solution file
  //              else: update the history file 
  bool          failed = false;
  std::ofstream fout;
  
  if ( is_sol )
    fout.open ( file_name.c_str() );
  else
    fout.open ( file_name.c_str() , std::ios::app );

  if ( !fout.fail() ) {

    fout.setf      ( std::ios::fixed             );
    fout.precision ( NOMAD::DISPLAY_PRECISION_BB );

    // solution display:
    if ( is_sol ) {
      if ( _p.get_bb_input_include_seed() )
	fout << _p.get_seed() << std::endl;
      if ( _p.get_bb_input_include_tag() )
	fout << x.get_tag() << std::endl;
      x.Point::display ( fout , "\n" , -1 , -1 );
      fout << std::endl;
    }

    // history display:
    else {
      x.Point::display ( fout , " " , -1 , -1 );
      fout << " ";
      x.get_bb_outputs().Point::display ( fout , " " , -1 , -1 );
      fout << std::endl;
    }

    if ( fout.fail() )
      failed = true;
  }
  else
    failed = true;

  fout.close();

  if ( failed && _p.out().get_gen_dd() != NOMAD::NO_DISPLAY )
    _p.out() << std::endl
	     << "Warning (" << "Evaluator_Control.cpp" << ", " << __LINE__
	     << "): could not "
	     << ( is_sol ? "save the current solution" :
		  "update the history" )
	     << " in \'"
	     << file_name << "\'" << std::endl << std::endl;
}

/*---------------------------------------------------------*/
/*             display evaluation result (private)         */
/*---------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_eval_result
( const NOMAD::Eval_Point & x                ,
  NOMAD::dd_type            display_degree   ,
  NOMAD::search_type        search           ,
  NOMAD::success_type       one_eval_success ,
  NOMAD::success_type       success            ) const
{
  const NOMAD::Display & out = _p.out();

  // surrogate evaluation:
  if ( x.get_eval_type() == NOMAD::SGTE ) {
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      out << std::endl << "point #" << x.get_tag() << " sgte eval: ";
      if ( x.is_eval_ok() ) {
	out << "h=";
	if ( x.get_h().is_defined() )
	  out << x.get_h();
	else
	  out << "inf (extr. barrier)";
	out << " f=" << x.get_f();
      }
      else
	out << "failed";
      out << std::endl;
    }
    if ( !_p.get_opt_only_sgte() )
      return;
  }

  // update the history file:
  // (contains surrogate evaluations if opt_only_sgte==true)
  // (history file is disabled during VNS search)
  const std::string & his_file = _p.get_history_file();
  if ( !his_file.empty() )
    write_sol_or_his_file ( _p.get_problem_dir() + his_file , x , false );

  // success displays:
  if ( one_eval_success != NOMAD::UNSUCCESSFUL &&
       one_eval_success >= success ) {

    // save the solution file:
    write_solution_file ( x );

    bool ds_ok = one_eval_success == NOMAD::FULL_SUCCESS &&
                 x.is_feasible ( _p.get_h_min() );

    // normal display:
    if ( display_degree == NOMAD::NORMAL_DISPLAY && ds_ok )
      display_stats ( false , out , _p.get_display_stats() , &x , NULL );

    // detailed display:
    else if ( display_degree == NOMAD::FULL_DISPLAY )
      out << std::endl << search << " " << one_eval_success
	  << " point " << x;

    // stats file:
    if ( ds_ok ) {
      const std::string & stats_file_name = _p.get_stats_file_name();
      if ( !stats_file_name.empty() )
	stats_file ( stats_file_name , &x , NULL );
    }
    
  }
  else if ( display_degree == NOMAD::FULL_DISPLAY ) {
    out << search << " " << one_eval_success
	<< " point #" << x.get_tag();
    if ( x.is_eval_ok() )
      out << " [ h=" << x.get_h()
	  << " f=" << x.get_f() << " ]" << std::endl;
    else
      out << ": evaluation failed" << std::endl;
  }
}

/*-------------------------------------------*/
/*        search a point in the cache        */
/*-------------------------------------------*/
/* . return true if the point is in cache    */
/* . private method                          */
/*-------------------------------------------*/
bool NOMAD::Evaluator_Control::cache_check
( const NOMAD::Eval_Point *& x              ,
  NOMAD::Barrier           & true_barrier   ,
  NOMAD::Barrier           & sgte_barrier   ,
  NOMAD::Pareto_Front      * pareto_front   ,
  bool                     & count_eval     ,
  const NOMAD::Double      & h_max          ,
  NOMAD::dd_type             display_degree   ) const
{
  NOMAD::eval_type          x_eval_type = x->get_eval_type();
  const NOMAD::Eval_Point * cache_x     = NULL;

  // first cache check:
  if ( x->is_in_cache() )
    cache_x = x;

  // second cache check:
  else
    cache_x = ( ( x->get_eval_type() == NOMAD::TRUTH ) ?
		_cache : _sgte_cache )->find ( *x );

  // cache hit: we transfer some data from x to cache_x:
  if ( cache_x ) {
     
    if ( x_eval_type != cache_x->get_eval_type() )
      throw NOMAD::Exception ( "Evaluator_Control.cpp" , __LINE__ ,
      "Evaluator_Control::cache_check(): eval and cache pts have different eval_type" );
 
    if ( cache_x->is_eval_ok() ) {

      NOMAD::Eval_Point * modifiable_cache_x
	= &NOMAD::Cache::get_modifiable_point ( *cache_x );

      // if wrong number of outputs, we reset cache_x._bb_outputs:
      {
	int m = _p.get_bb_nb_outputs();
	if ( cache_x->get_bb_outputs().size() != m )
	  modifiable_cache_x->set_bb_output ( NOMAD::Point ( m ) );
      }
	
      modifiable_cache_x->set_signature          ( x->get_signature         () );
      modifiable_cache_x->set_direction          ( x->get_direction         () );
      modifiable_cache_x->set_mesh_index         ( x->get_mesh_index        () );
      modifiable_cache_x->set_poll_center_type   ( x->get_poll_center_type  () );
      modifiable_cache_x->set_user_eval_priority ( x->get_user_eval_priority() );
      modifiable_cache_x->set_rand_eval_priority ( x->get_rand_eval_priority() );

      // set_f, set_h, and set_EB_ok:
      _ev->compute_f ( *modifiable_cache_x );
      _ev->compute_h ( *modifiable_cache_x );
    }
  }
  
  // point in cache but evaluation is to be made again:
  if ( cache_x && cache_x->is_eval_ok() && 
       ( !cache_x->get_f().is_defined() ||
	 ( cache_x->is_EB_ok()                      &&
	   !cache_x->get_bb_outputs().is_complete() &&
	   cache_x->get_h().is_defined()            &&
	   cache_x->get_h() < h_max                    ) ) ) {
    x       = cache_x;
    cache_x = NULL;
  }

  // point in cache:
  if ( cache_x ) {

    _stats.add_cache_hit();

    // displays:
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      const NOMAD::Display & out = _p.out();
      if ( cache_x->get_eval_type() == NOMAD::SGTE )
	out << "surrogate ";
      out << "cache hit: #" << x->get_tag()
	  << " --> #" << cache_x->get_tag() << std::endl;
    }

    // we process the Eval_Point taken in cache:
    process_eval_point ( *cache_x ,
			 ( cache_x->get_eval_type() == NOMAD::TRUTH ) ?
			 true_barrier : sgte_barrier ,
			 pareto_front );

    // count the (simulated) bb eval ?
    int index_cnt_eval = _p.get_index_cnt_eval();
    if ( index_cnt_eval >= 0 && cache_x->get_bb_outputs()[index_cnt_eval] == 0.0 )
      count_eval = false;

    x = cache_x;

    return true;
  }

  return false;
}

/*----------------------------------------------------*/
/*                 eval a point (private)             */
/*----------------------------------------------------*/
void NOMAD::Evaluator_Control::eval_point ( const NOMAD::Eval_Point * x            ,
					    NOMAD::Barrier          & true_barrier ,
					    NOMAD::Barrier          & sgte_barrier ,
					    NOMAD::Pareto_Front     * pareto_front ,
					    bool                    & count_eval   ,
					    bool                    & stop         ,
					    NOMAD::stop_type        & stop_reason  ,
					    const NOMAD::Double     & h_max          )
{
  int max_bb_eval   = _p.get_max_bb_eval();
  int max_sgte_eval = _p.get_max_sgte_eval();

  // blackbox or surrogate evaluations are allowed:
  if ( ( x->get_eval_type() == NOMAD::TRUTH && max_bb_eval   != 0 ) ||
       ( x->get_eval_type() == NOMAD::SGTE  && max_sgte_eval != 0 )    ) {
    
    NOMAD::Eval_Point * eval_x = &NOMAD::Cache::get_modifiable_point ( *x );

    // get the signature:
    NOMAD::Signature * signature = x->get_signature();
    if ( !signature )
      throw NOMAD::Exception ( "Evaluator_Control.cpp" , __LINE__ ,
	    "Evaluator_Control::eval_point(): the point has no signature" );

    // evaluation of the point:
    // ------------------------
    bool eval_ok = true;

    {
      // 1. scaling:
      bool do_scaling = signature->get_scaling().is_defined();
      if ( do_scaling )
	eval_x->scale();

      // 2. evaluation:
      try {
	eval_ok = _ev->eval_x ( *eval_x , h_max , count_eval );
      }
      catch ( ... ) {
	eval_ok = false;
      }
    
      // 3. unscaling:
      if ( do_scaling )
	eval_x->unscale();
    }

    if ( eval_ok ) {
      
      eval_x->set_eval_status ( NOMAD::EVAL_OK );
      
      // set_f, set_h and set_EB_ok:
      _ev->compute_f ( *eval_x );
      _ev->compute_h ( *eval_x );
      
      // we process the evaluated point:
      process_eval_point ( *eval_x ,
 			   ( eval_x->get_eval_type() == NOMAD::TRUTH ) ?
 			   true_barrier : sgte_barrier ,
 			   pareto_front );
    }
    else {
      eval_x->set_eval_status ( NOMAD::EVAL_FAIL );
      _stats.add_failed_eval();
    }

    // insertion in cache even if is_eval_ok == false:
    if ( !x->is_in_cache() )
      ( ( x->get_eval_type() == NOMAD::SGTE ) ?
	_sgte_cache : _cache)->insert ( *x );
    
    // we count the bb evaluation:
    if ( count_eval ) {
      if ( x->get_eval_type() == NOMAD::SGTE )
	_stats.add_sgte_eval();
      else
	_stats.add_bb_eval();
    }
    
    // we count the output stats (STAT_SUM and STAT_AVG):
    if ( _p.check_stat_sum() || _p.check_stat_avg() ) {
      
      count_output_stats ( *x );
      
      // check STAT_SUM_TARGET:
      NOMAD::Double sum_target = _p.get_stat_sum_target();
      if ( sum_target.is_defined() ) {
	NOMAD::Double sum = _stats.get_stat_sum();
	if ( !stop && sum.is_defined() && sum >= sum_target ) {
	  stop        = true;
	  stop_reason = NOMAD::STAT_SUM_TARGET_REACHED;
	}
      }
    }
  }
  
  // check the number of blackbox evaluations:
  if ( !stop ) {
    if ( max_bb_eval > 0 && _stats.get_bb_eval() >= max_bb_eval ) {
      stop        = true;
      stop_reason = NOMAD::MAX_BB_EVAL_REACHED;
    }
    if ( max_sgte_eval > 0 && _stats.get_sgte_eval() >= max_sgte_eval ) {
      stop        = true;
      stop_reason = NOMAD::MAX_SGTE_EVAL_REACHED;
    }
  } 
}

/*-------------------------------------------*/
/*      check stopping criteria (private)    */
/*-------------------------------------------*/
void NOMAD::Evaluator_Control::check_stopping_criteria
( NOMAD::search_type      search ,
  bool                count_eval ,
  const NOMAD::Eval_Point    * x ,
  bool                    & stop ,
  NOMAD::stop_type & stop_reason   ) const
{
  // check the time:
  if ( !stop                 &&
       _p.get_max_time() > 0 &&
       _stats.get_real_time() >= _p.get_max_time() ) {
    stop        = true;
    stop_reason = NOMAD::MAX_TIME_REACHED;
  }

  // count an evaluation or a simulated blackbox evaluation:
  if ( x->get_eval_type() == NOMAD::TRUTH ) {
    _stats.add_eval();
    if ( count_eval && !x->get_current_run() )
      _stats.add_sim_bb_eval();
  }

  // check the stopping condition MAX_EVAL:
  if ( !stop                 &&
       _p.get_max_eval() > 0 &&
       _stats.get_eval() >= _p.get_max_eval() ) {
    stop        = true;
    stop_reason = NOMAD::MAX_EVAL_REACHED;
  }
  
  // check the stopping condition MAX_SIM_BB_EVAL:
  if ( !stop                         &&
       _p.get_max_sim_bb_eval() >  0 &&
       _stats.get_sim_bb_eval() >= _p.get_max_sim_bb_eval() ) {
    stop        = true;
    stop_reason = NOMAD::MAX_SIM_BB_EVAL_REACHED;
  }
  
  // check the stopping conditions F_TARGET and FEAS_REACHED
  // (for phase one: the evaluations must stop if all EB
  //  constraints are satisfied, but some PB constraints can
  //  be violated)
  if ( !stop           &&
       x->is_eval_ok() &&
       ( _p.get_opt_only_sgte() ||
	 x->get_eval_type() == NOMAD::TRUTH ) ) {
    
    bool feasible = x->is_feasible ( _p.get_h_min() );
    
    // check FEAS_REACHED:
    if ( feasible && _p.get_stop_if_feasible() ) {	  
      stop        = true;
      stop_reason = NOMAD::FEAS_REACHED;
    }
    
    // check F_TARGET:
    {
      const NOMAD::Point           & f_target       = _p.get_f_target();
      const std::list<int>         & index_obj      = _p.get_index_obj();
      std::list<int>::const_iterator index_obj_end  = index_obj.end();
      bool                           check_f_target = f_target.is_defined();
      int                            nb_to_check    = (check_f_target) ?
	                                              f_target.nb_defined() : 0;
      
      if ( check_f_target && ( feasible || search == NOMAD::LH_SEARCH_P1 ) ) {
	const NOMAD::Point & bbo = x->get_bb_outputs();
	bool                 chk = true;
	int                  k   = 0;
	int                  cnt = 0;
	for ( std::list<int>::const_iterator it = index_obj.begin();
	      it != index_obj_end ; ++it , ++k ) {
	  if ( bbo[*it].is_defined() && f_target[k].is_defined() ) {
	    if ( f_target[k] < bbo[*it] ) {
	      chk = false;
	      break;
	    }
	    cnt++;
	  }
	}
	
	if ( chk && cnt == nb_to_check ) {
	  stop        = true;
	  stop_reason = NOMAD::F_TARGET_REACHED;
	}
      }
    }
  }
}

/*-------------------------------------------------------*/
/*  receive an evaluation result from a slave (private)  */
/*-------------------------------------------------------*/
#ifdef USE_MPI
void NOMAD::Evaluator_Control::receive_eval_result
( NOMAD::search_type    search       ,
  NOMAD::Eval_Point   * x            ,
  NOMAD::Barrier      & true_barrier ,
  NOMAD::Barrier      & sgte_barrier ,
  NOMAD::Pareto_Front * pareto_front ,
  int                   slave_rank   ,
  bool                & stop         ,
  NOMAD::stop_type    & stop_reason    )
{
  bool eval_ok , count_eval;

  // receive the evaluation result:
  _slave->receive_eval_result ( slave_rank , x , eval_ok , count_eval );

  // process the evaluation:
  if ( eval_ok ) {
      
    // set_f, set_h and set_EB_ok:
    _ev->compute_f ( *x );
    _ev->compute_h ( *x );
      
    // we process the evaluated point:
    process_eval_point ( *x                                    ,
			 (x->get_eval_type()==NOMAD::TRUTH) ?
			 true_barrier : sgte_barrier           ,
			 pareto_front                            );
  }
  else
    _stats.add_failed_eval();
  
  // insertion in cache even if !eval_ok:
  if ( !x->is_in_cache() )
    ( ( x->get_eval_type() == NOMAD::SGTE ) ?
      _sgte_cache : _cache)->insert ( *x );
    
  // we count the bb evaluation:
  if ( count_eval ) {
    if ( x->get_eval_type() == NOMAD::SGTE )
      _stats.add_sgte_eval();
    else
      _stats.add_bb_eval();
  }
    
  // we count the output stats (STAT_SUM and STAT_AVG):
  if ( _p.check_stat_sum() || _p.check_stat_avg() ) {
      
    count_output_stats ( *x );
      
    // check STAT_SUM_TARGET:
    NOMAD::Double sum_target = _p.get_stat_sum_target();
    if ( sum_target.is_defined() ) {
      NOMAD::Double sum = _stats.get_stat_sum();
      if ( !stop && sum.is_defined() && sum >= sum_target ) {
	stop        = true;
	stop_reason = NOMAD::STAT_SUM_TARGET_REACHED;
      }
    }
  }

  // check stopping criteria:
  if ( !stop ) {

    int max_bb_eval   = _p.get_max_bb_eval();
    int max_sgte_eval = _p.get_max_sgte_eval();
    
    if ( max_bb_eval > 0 && _stats.get_bb_eval() >= max_bb_eval ) {
      stop        = true;
      stop_reason = NOMAD::MAX_BB_EVAL_REACHED;
    }
    if ( max_sgte_eval > 0 && _stats.get_sgte_eval() >= max_sgte_eval ) {
      stop        = true;
      stop_reason = NOMAD::MAX_SGTE_EVAL_REACHED;
    }
  }

  check_stopping_criteria ( search , count_eval , x , stop , stop_reason );
}
#endif

/*----------------------------------------------------*/
/*          wait for evaluations in progress          */
/*----------------------------------------------------*/
#ifdef USE_MPI
void NOMAD::Evaluator_Control::wait_for_evaluations
( NOMAD::search_type                     search         ,
  NOMAD::Barrier                       & true_barrier   ,
  NOMAD::Barrier                       & sgte_barrier   ,
  NOMAD::Pareto_Front                  * pareto_front   ,
  bool                                 & stop           ,
  NOMAD::stop_type                     & stop_reason    ,
  NOMAD::success_type                  & success        ,
  std::list<const NOMAD::Eval_Point *> & evaluated_pts    )
{ 
  if ( _nb_in_progress == 0 )
    return;

  // display degree:
  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_display_degree ( search );

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl
	<< NOMAD::open_block ( "wait for evaluations" );

  NOMAD::Barrier     & barrier = ( _p.get_opt_only_sgte() ) ?
                                 sgte_barrier : true_barrier;
  char                 signal;
  int                  source;
  NOMAD::Eval_Point  * eval_x;
  NOMAD::success_type  one_eval_success;

  while ( _nb_in_progress > 0 ) {
  
    source = NOMAD::Slave::receive_signal ( signal );  
    eval_x = _eval_in_progress[source];

    if ( eval_x ) {

      if ( display_degree == NOMAD::FULL_DISPLAY )
	out << std::endl << "receive eval point #" << eval_x->get_tag()
	    << " from slave " << source << std::endl << std::endl;

      receive_eval_result ( search       ,
			    eval_x       ,
			    true_barrier ,
			    sgte_barrier ,
			    pareto_front ,
			    source       ,
			    stop         ,
			    stop_reason    );
	
      // list of processed points:
      if ( eval_x->is_in_cache() )
	evaluated_pts.push_back ( eval_x );
	
      // success:
      one_eval_success = barrier.get_one_eval_succ();
      success          = barrier.get_success();
      
      // asynchronous success count:
      if ( success == NOMAD::FULL_SUCCESS &&
	   _elop_tag != _slaves_elop_tags[source] )
	_stats.add_asynchronous_success();

      // displays:
      display_eval_result ( *eval_x          ,
			    display_degree   ,
			    search           ,
			    one_eval_success ,
			    success            );
      
      if ( !_eval_in_progress[source]->is_in_cache() )
	delete _eval_in_progress[source];
      _eval_in_progress[source] = NULL;
      _slaves_elop_tags[source] = -1;
      --_nb_in_progress;
      
      // force quit (by pressing ctrl-c):
      if ( !stop && NOMAD::Evaluator_Control::_force_quit ) {
	stop        = true;
	stop_reason = NOMAD::CTRL_C;
	break;
      }
      
      if ( stop && ( stop_reason==NOMAD::ERROR ||
		     stop_reason==NOMAD::UNKNOWN_STOP_REASON ) )
	break;
    }
    else
      NOMAD::Slave::send_signal ( NOMAD::WAIT_SIGNAL , source );
  }

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out.close_block();
}
#endif

/*----------------------------------------------------------------*/
/*  check if the evaluation at this point is already in progress  */
/*  (private)                                                     */
/*----------------------------------------------------------------*/
#ifdef USE_MPI
bool NOMAD::Evaluator_Control::already_in_progress
( const NOMAD::Eval_Point & x ) const
{
  if ( _eval_in_progress ) {

    int x_tag = x.get_tag();
    int np    = NOMAD::Slave::get_nb_processes();

    for ( int i = 0 ; i < np ; ++i )
      if ( _eval_in_progress[i] &&
	   ( _eval_in_progress[i]->get_tag() == x_tag ||
	     _eval_in_progress[i]->Point::operator == ( x ) ) )
	return true;
  }
  return false;
}
#endif

/*----------------------------------------------------------------*/
/*     eval_list_of_points, private version (parallel version)    */
/*----------------------------------------------------------------*/
#ifdef USE_MPI
void NOMAD::Evaluator_Control::private_eval_list_of_points
( NOMAD::search_type              search         ,   // IN     : search type
  NOMAD::Barrier                & true_barrier   ,   // IN/OUT : the barrier
  NOMAD::Barrier                & sgte_barrier   ,   // IN/OUT : the surrogate barrier
  NOMAD::Pareto_Front           * pareto_front   ,   // IN/OUT : the Pareto front
                                                     //          (can be NULL)
  bool                          & stop           ,   // IN/OUT : stopping criterion
  NOMAD::stop_type              & stop_reason    ,   // OUT    : stopping reason
  const NOMAD::Eval_Point      *& new_feas_inc   ,   // OUT    : new feasible incumbent
  const NOMAD::Eval_Point      *& new_infeas_inc ,   // OUT    : new infeas. incumbent
  NOMAD::success_type           & success        ,   // OUT    : type of success
  std::list<const NOMAD::Eval_Point *>
                                & evaluated_pts    ) // OUT    : list of processed pts
{
  if ( stop || _eval_lop.empty() ) {
    stop_reason = NOMAD::UNKNOWN_STOP_REASON;
    ++_elop_tag;
    return;
  }

  evaluated_pts.clear();

  // initial display:
  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_display_degree ( search );

  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    std::ostringstream msg;
    msg << "list of points evaluation (" << search << ")";
    out << std::endl << NOMAD::open_block ( msg.str() );
  }

  // call the Evaluator (virtual) preprocessing of a list of points:
  _ev->list_of_points_preprocessing ( _eval_lop );

  const NOMAD::Eval_Point * old_feasible_incumbent   = NULL;
  const NOMAD::Eval_Point * old_infeasible_incumbent = NULL;

  // active barrier:
  NOMAD::Barrier & barrier = ( _p.get_opt_only_sgte() ) ?
    sgte_barrier : true_barrier;

  old_feasible_incumbent   = barrier.get_best_feasible();
  old_infeasible_incumbent = barrier.get_best_infeasible();
  
  NOMAD::Double f0;
  if ( _p.get_opportunistic_min_f_imprvmt().is_defined() &&
       old_feasible_incumbent )
    f0 = old_feasible_incumbent->get_f();

  new_feas_inc   = NULL;
  new_infeas_inc = NULL;
  stop           = false;
  success        = NOMAD::UNSUCCESSFUL;
  stop_reason    = NOMAD::NO_STOP;
  
  const NOMAD::Eval_Point  * x;
  NOMAD::check_failed_type   check_failed_reason;
  bool                       count_eval;
  std::vector<const NOMAD::Eval_Point *>
                             to_be_evaluated;
  NOMAD::success_type        one_eval_success;
  bool                       one_for_luck = false;
  bool                       opp_stop     = false;
  int                        init_nb_eval = _stats.get_eval();
  int                        nb_success   = 0;
  int                        k            = 0;
  int                        nb_points    = static_cast<int> ( _eval_lop.size() );
  int                        max_bb_eval  = _p.get_max_bb_eval();

  // loop #1: search in cache:
  // -------------------------
  std::set<NOMAD::Priority_Eval_Point>::iterator
    it  = _eval_lop.begin() ,
    end = _eval_lop.end();
  while ( !stop && !opp_stop && it != end ) {

    x = it->get_point();

    x->set_current_run ( true );

    // displays:
    if ( display_degree == NOMAD::FULL_DISPLAY ) {

      // open the evaluation block:
      {
	std::ostringstream oss;
	if ( x->get_eval_type() == NOMAD::SGTE )
	  oss << "surrogate ";
	oss << "evaluation " << k+1 << "/" << nb_points;
	out << std::endl << NOMAD::open_block ( oss.str() );
      }

      out << std::endl << "point #" << x->get_tag() << " ( ";
      x->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
      out << " )" << std::endl;
      if ( x->get_direction() )
 	out << "direction  : " << *x->get_direction()  << std::endl;
      if ( x->get_mesh_index() )
 	out << "mesh index : " << *x->get_mesh_index() << std::endl;
    }

    // check if the evaluation at this point is already in progress:
    if ( !already_in_progress ( *x ) ) {

      // current point check (# of bb outputs, bounds, integer values, fixed-vars):
      if ( x->check ( _p.get_bb_nb_outputs() , check_failed_reason ) ) {

	count_eval = true;

	// point in cache:
	if ( cache_check ( x                   ,
			   true_barrier        ,
			   sgte_barrier        ,
			   pareto_front        ,
			   count_eval          ,
			   barrier.get_h_max() ,
			   display_degree        ) ) {
   
	  // list of processed points:
	  evaluated_pts.push_back ( x );

	  // check stopping criteria:
	  check_stopping_criteria ( search , count_eval , x , stop , stop_reason );

	  // success:
	  one_eval_success = barrier.get_one_eval_succ();
	  success          = barrier.get_success();
	  
	  // displays:
	  display_eval_result ( *x               ,
				display_degree   ,
				search           ,
				one_eval_success ,
				success            );
	  
	  // stop the evaluations (opportunistic strategy) ?
	  if ( stop_evaluations ( *x               ,
				  search           ,
				  k                ,
				  nb_points        ,
				  stop             ,
				  display_degree   ,
				  one_eval_success ,
				  success          ,
				  init_nb_eval     ,
				  f0               ,
				  barrier          ,
				  nb_success       ,
				  one_for_luck       ) ) {
	    _stats.add_interrupted_eval();
	    opp_stop = true; // will break loop #1
	  }

	  // close the evaluation block:
	  if ( display_degree == NOMAD::FULL_DISPLAY )
	    out.close_block();
	}
	
	// point not in cache (the point is saved for loop #2):
	else {
	  
	  // blackbox or surrogate evaluations are allowed:
	  if ( ( x->get_eval_type() == NOMAD::TRUTH && max_bb_eval != 0 ) ||
	       ( x->get_eval_type() == NOMAD::SGTE  && _p.get_max_sgte_eval() != 0 ) )
	    to_be_evaluated.push_back ( x );

	  // close the evaluation block:
	  if ( display_degree == NOMAD::FULL_DISPLAY )
	    out.close_block();
	}
      }
      
      // points[k]->check() failed (close the evaluation block):
      else if ( display_degree == NOMAD::FULL_DISPLAY ) {
	std::ostringstream oss;
	oss << "check failed (" << check_failed_reason << ")";
	out.close_block ( oss.str() );
      }      
    }
    
    // evaluation already in progress (close the evaluation block):
    else if ( display_degree == NOMAD::FULL_DISPLAY ) {
      std::ostringstream oss;
      oss << "evaluation of point #" << x->get_tag()
	  << " already in progress";
      out.close_block ( oss.str() );
    }
    
    ++it;
    ++k;

    // force quit (by pressing ctrl-c):
    if ( !stop && NOMAD::Evaluator_Control::_force_quit ) {
      stop        = true;
      stop_reason = NOMAD::CTRL_C;
    }

  }  // end of loop #1
     // --------------

  // loop #2: evaluations:
  // ---------------------
  int                 nb_to_evaluate = static_cast<int> ( to_be_evaluated.size() );
  int                 nb_evaluated   = 0;
  int                 cur            = 0;
  int                 source;
  char                signal;
  NOMAD::Eval_Point * eval_x;

  while ( !stop && !opp_stop && nb_evaluated < nb_to_evaluate ) {

    source = NOMAD::Slave::receive_signal ( signal );

    // 2.1: send the RESULT signal, receive and process the evaluation result:
    // -----------------------------------------------------------------------
    eval_x = _eval_in_progress[source];
    if ( eval_x ) {
      
      if ( display_degree == NOMAD::FULL_DISPLAY )
	out << std::endl << "receive eval point #" << eval_x->get_tag()
	    << " from slave " << source << std::endl << std::endl;

      receive_eval_result ( search       ,
			    eval_x       ,
			    true_barrier ,
			    sgte_barrier ,
			    pareto_front ,
			    source       ,
			    stop         ,
			    stop_reason    );

      // list of processed points:
      if ( eval_x->is_in_cache() )
	evaluated_pts.push_back ( eval_x );
            
      // success:
      one_eval_success = barrier.get_one_eval_succ();
      success          = barrier.get_success();
      
      // asynchronous success count:
      if ( success == NOMAD::FULL_SUCCESS &&
	   _elop_tag != _slaves_elop_tags[source] )
	_stats.add_asynchronous_success();

      // displays:
      display_eval_result ( *eval_x          ,
			    display_degree   ,
			    search           ,
			    one_eval_success ,
			    success            );
      
      // stop the evaluations (opportunistic strategy) ?
      if ( stop_evaluations ( *eval_x          ,
			      search           ,
			      nb_evaluated     ,
			      nb_to_evaluate   ,
			      stop             ,
			      display_degree   ,
			      one_eval_success ,
			      success          ,
			      init_nb_eval     ,
			      f0               ,
			      barrier          ,
			      nb_success       ,
			      one_for_luck       ) ) {
	_stats.add_interrupted_eval();
	opp_stop = true; // will break loop #2
      }
      
      _eval_in_progress[source] = NULL;
      _slaves_elop_tags[source] = -1;
      --_nb_in_progress;
      ++nb_evaluated;
    }
    
    // 2.2: send the EVAL signal and launch a new evaluation:
    // ------------------------------------------------------
    else {
      
      // do not launch a new evaluation if...

      // there is no more points to be evaluated:
      if ( cur == nb_to_evaluate )
	NOMAD::Slave::send_signal ( NOMAD::WAIT_SIGNAL , source );

      // or if bbe+_nb_in_progress >= max_bb_eval:
      else if ( to_be_evaluated[cur]->get_eval_type() == NOMAD::TRUTH &&
		max_bb_eval > 0 &&
		_stats.get_bb_eval() + _nb_in_progress >= max_bb_eval    ) {
	stop        = true;
	stop_reason = NOMAD::MAX_BB_EVAL_REACHED;
	NOMAD::Slave::send_signal ( NOMAD::WAIT_SIGNAL , source );
      }

      else {
	
	// get the signature:
	NOMAD::Signature * signature = to_be_evaluated[cur]->get_signature();

	// there is no signature (error):
	if ( !signature ) {
	  stop        = true;
	  stop_reason = NOMAD::ERROR;
	  if ( display_degree != NOMAD::NO_DISPLAY )
	    out << std::endl
		<< "Error in Evaluator_Control::private_eval_list_of_points():"
		<< " the point #" << to_be_evaluated[cur]->get_tag()
		<< " has no signature" << std::endl << std::endl;
	  NOMAD::Slave::send_signal ( NOMAD::WAIT_SIGNAL , source );
	}
	
	else {

	  NOMAD::Slave::send_signal ( NOMAD::EVAL_SIGNAL , source );

	  eval_x = &NOMAD::Cache::get_modifiable_point ( *to_be_evaluated[cur++] );

	  if ( display_degree == NOMAD::FULL_DISPLAY )
	    out << std::endl
		<< "send eval point #" << eval_x->get_tag()
		<< " to slave " << source << std::endl;
      
	  // 1. scaling:
	  bool do_scaling = signature->get_scaling().is_defined();
	  if ( do_scaling )
	    eval_x->scale();

	  // 2. send the point:
	  _slave->send_eval_point ( eval_x , source , barrier.get_h_max() );
      
	  // 3. unscaling:
	  if ( do_scaling )
	    eval_x->unscale();

	  eval_x->set_eval_status ( NOMAD::EVAL_IN_PROGRESS );
	
	  _eval_in_progress[source] = eval_x;
	  _slaves_elop_tags[source] = _elop_tag;
	  ++_nb_in_progress;
	}
      }
    }
    
    // force quit (by pressing ctrl-c):
    if ( !stop && NOMAD::Evaluator_Control::_force_quit ) {
      stop        = true;
      stop_reason = NOMAD::CTRL_C;
    }

  }  // end of loop #2
     // --------------

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl
	<< "number of evaluations in progress: " << _nb_in_progress
	<< std::endl << std::endl;

  // the algorithm is not asynchronous: we have to wait for
  // all the evaluations in progress:
  if ( !_p.get_asynchronous() )
    
    wait_for_evaluations ( search        ,
 			   true_barrier  ,
 			   sgte_barrier  ,
 			   pareto_front  ,
 			   stop          ,
			   stop_reason   ,
			   success       ,
			   evaluated_pts   );
    
  // barriers update:
  if ( !stop ) {
    true_barrier.update_and_reset_success();
    sgte_barrier.update_and_reset_success();
  }
  
  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << NOMAD::close_block ( "end of evaluations" ) << std::endl;
  
  // incumbents update:
  const NOMAD::Eval_Point * bf = barrier.get_best_feasible  ();
  const NOMAD::Eval_Point * bi = barrier.get_best_infeasible();
  if ( bf && bf != old_feasible_incumbent )
    new_feas_inc = bf;
  if ( bi && bi != old_infeasible_incumbent )
    new_infeas_inc = bi;
  
  // the list of eval. points is deleted (only points in the cache are kept):
  clear_eval_lop();
  
  // update the unique eval_lop() tag:
  ++_elop_tag;
  
} // end of eval_lop() parallel version

/*----------------------------------------------------------------*/
/*       eval_list_of_points, private version (scalar version)    */
/*----------------------------------------------------------------*/
#else
void NOMAD::Evaluator_Control::private_eval_list_of_points
( NOMAD::search_type              search         ,   // IN     : search type
  NOMAD::Barrier                & true_barrier   ,   // IN/OUT : the barrier
  NOMAD::Barrier                & sgte_barrier   ,   // IN/OUT : the surrogate barrier
  NOMAD::Pareto_Front           * pareto_front   ,   // IN/OUT : the Pareto front
                                                     //          (can be NULL)
  bool                          & stop           ,   // IN/OUT : stopping criterion
  NOMAD::stop_type              & stop_reason    ,   // OUT    : stopping reason
  const NOMAD::Eval_Point      *& new_feas_inc   ,   // OUT    : new feasible incumbent
  const NOMAD::Eval_Point      *& new_infeas_inc ,   // OUT    : new infeas. incumbent
  NOMAD::success_type           & success        ,   // OUT    : type of success
  std::list<const NOMAD::Eval_Point *>
                                & evaluated_pts    ) // OUT    : list of processed pts
{
  if ( stop || _eval_lop.empty() ) {
    stop_reason = NOMAD::UNKNOWN_STOP_REASON;
    return;
  }

  evaluated_pts.clear();

  // initial display:
  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_display_degree ( search );

  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    std::ostringstream oss;
    oss << "list of points evaluation (" << search << ")";
    out << std::endl << NOMAD::open_block ( oss.str() );
  }

  // call the Evaluator (virtual) preprocessing of a list of points:
  _ev->list_of_points_preprocessing ( _eval_lop );

  const NOMAD::Eval_Point * old_feasible_incumbent   = NULL;
  const NOMAD::Eval_Point * old_infeasible_incumbent = NULL;

  // active barrier:
  NOMAD::Barrier & barrier = ( _p.get_opt_only_sgte() ) ?
    sgte_barrier : true_barrier;

  old_feasible_incumbent   = barrier.get_best_feasible();
  old_infeasible_incumbent = barrier.get_best_infeasible();
  
  NOMAD::Double f0;
  if ( _p.get_opportunistic_min_f_imprvmt().is_defined() &&
       old_feasible_incumbent )
    f0 = old_feasible_incumbent->get_f();

  new_feas_inc   = NULL;
  new_infeas_inc = NULL;
  stop           = false;
  success        = NOMAD::UNSUCCESSFUL;
  stop_reason    = NOMAD::NO_STOP;
  
  const NOMAD::Eval_Point * x;
  NOMAD::check_failed_type  check_failed_reason;
  bool                      count_eval;
  bool                      one_for_luck = false;
  bool                      stop_evals   = false;
  int                       init_nb_eval = _stats.get_eval();
  int                       nb_success   = 0;
  int                       k            = 0;
  int                       nb_points    = static_cast<int>(_eval_lop.size());

  // main loop (on the list of points):
  // ----------------------------------
  std::set<NOMAD::Priority_Eval_Point>::iterator
    it  = _eval_lop.begin() ,
    end = _eval_lop.end();

  while ( !stop_evals && !stop && it != end ) {

    x = it->get_point();
    
    x->set_current_run ( true );

    // displays:
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      {
	// open the evaluation block:
	std::ostringstream oss;
	if ( x->get_eval_type() == NOMAD::SGTE )
	  oss << "surrogate ";
	oss << "evaluation " << k+1 << "/" << nb_points;
	out << std::endl << NOMAD::open_block ( oss.str() );
      }

      out << std::endl << "point #" << x->get_tag() << " ( ";
      x->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
      out << " )" << std::endl;
      if ( x->get_direction() )
	out << "direction  : " << *x->get_direction()  << std::endl;
      if ( x->get_mesh_index() )
	out << "mesh index : " << *x->get_mesh_index() << std::endl;
      out << std::endl;
    }
    
    // current point check (# of bb outputs, bounds, integer values, fixed-vars):
    if ( x->check ( _p.get_bb_nb_outputs() , check_failed_reason ) ) {

      count_eval = true;

      // search in cache or eval the point:
      if ( !cache_check ( x                   ,
			  true_barrier        ,
			  sgte_barrier        ,
			  pareto_front        ,
			  count_eval          ,
			  barrier.get_h_max() ,
			  display_degree        ) )
	eval_point ( x                   ,
		     true_barrier        ,
		     sgte_barrier        ,
		     pareto_front        ,
		     count_eval          ,
		     stop                ,
		     stop_reason         ,
		     barrier.get_h_max()   );

      // list of processed points:
      if ( x->is_in_cache() )
	evaluated_pts.push_back ( x );

      // check stopping criteria:
      check_stopping_criteria ( search , count_eval , x , stop , stop_reason );

      // success:
      NOMAD::success_type one_eval_success = barrier.get_one_eval_succ();
      success                              = barrier.get_success();

      // displays:
      display_eval_result ( *x, display_degree, search, one_eval_success, success );

      // stop the evaluations (opportunistic strategy) ?
      if ( stop_evaluations ( *x               ,
			      search           ,
			      k                ,
			      nb_points        ,
			      stop             ,
			      display_degree   ,
			      one_eval_success ,
			      success          ,
			      init_nb_eval     ,
			      f0               ,
			      barrier          ,
			      nb_success       ,
			      one_for_luck       ) ) {
	_stats.add_interrupted_eval();
	stop_evals = true;
      }
    }

    // points[k]->check() failed:
    else if ( display_degree == NOMAD::FULL_DISPLAY )
      out << "check failed (" << check_failed_reason << ")" << std::endl;

    // close the evaluation block:
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << NOMAD::close_block();

    ++it;
    ++k;

    // force quit (by pressing ctrl-c):
    if ( !stop && NOMAD::Evaluator_Control::_force_quit ) {
      stop        = true;
      stop_reason = NOMAD::CTRL_C;
    }

  } // end of main loop
    // ----------------

  // barriers update:
  if ( !stop ) {
    true_barrier.update_and_reset_success();
    sgte_barrier.update_and_reset_success();
  }

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl << NOMAD::close_block ( "end of evaluations" )
	<< std::endl;

  // incumbents update:
  const NOMAD::Eval_Point * bf = barrier.get_best_feasible  ();
  const NOMAD::Eval_Point * bi = barrier.get_best_infeasible();
  if ( bf && bf != old_feasible_incumbent )
    new_feas_inc = bf;
  if ( bi && bi != old_infeasible_incumbent )
    new_infeas_inc = bi;
  
  // the list of eval. points is deleted (only points in the cache are kept):
  clear_eval_lop();

} // end of eval_lop() scalar version
#endif

/*-------------------------------------------------*/
/*  clear the list of evaluation points (private)  */
/*-------------------------------------------------*/
void NOMAD::Evaluator_Control::clear_eval_lop ( void )
{
  if ( _eval_lop.empty() )
    return;

  const NOMAD::Eval_Point * x;
  std::set<NOMAD::Priority_Eval_Point>::iterator
    it  = _eval_lop.begin() ,
    end = _eval_lop.end();

  while  ( it != end ) {
    x = it->get_point();

    if ( !x->is_in_cache() && x->get_eval_status() != NOMAD::EVAL_IN_PROGRESS )
      delete x;

    _eval_lop.erase ( it++ );
  }
}

/*----------------------------------------------------------------------------------*/
/*  evaluation of a list of points (public version that calls the private version)  */
/*----------------------------------------------------------------------------------*/
void NOMAD::Evaluator_Control::eval_list_of_points
( NOMAD::search_type              search             , // IN     : search type
  NOMAD::Barrier                & true_barrier       , // IN/OUT : truth barrier
  NOMAD::Barrier                & sgte_barrier       , // IN/OUT : surrogate barrier
  NOMAD::Pareto_Front           * pareto_front       , // IN/OUT : Pareto front
                                                       //          (can be NULL)
  bool                          & stop               , // IN/OUT : stopping criterion
  NOMAD::stop_type              & stop_reason        , // OUT    : stopping reason
  const NOMAD::Eval_Point      *& new_feas_inc       , // OUT    : new feas. incumbent
  const NOMAD::Eval_Point      *& new_infeas_inc     , // OUT    : new infeas. incumb.
  NOMAD::success_type           & success            , // OUT    : type of success
  std::list<const NOMAD::Eval_Point *>
                                * evaluated_pts   )    // OUT    : list of processed
                                                       //          pts (can be NULL)
{
  bool del_evaluated_pts = false;
  if ( !evaluated_pts ) {
    evaluated_pts     = new std::list<const NOMAD::Eval_Point *>;
    del_evaluated_pts = true;
  }

  bool sgte_eval_sort   = _p.get_sgte_eval_sort() && (_eval_lop.size() > 1);
  bool opt_only_sgte    = _p.get_opt_only_sgte();
  bool display_new_list = false;

  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_display_degree ( search );

  // reset the success type:
  true_barrier.reset_success();
  sgte_barrier.reset_success();

  // define all points as surrogates:
  if ( opt_only_sgte || sgte_eval_sort ) {
    for ( std::set<NOMAD::Priority_Eval_Point>::iterator it = _eval_lop.begin() ;
	  it != _eval_lop.end() ; ++it )
      NOMAD::Cache::get_modifiable_point(*it->get_point()).set_eval_type(NOMAD::SGTE);
  }

  // use the surrogates to sort the eval. points:
  if ( !opt_only_sgte && sgte_eval_sort ) {

    // evaluate the surrogate:
    private_eval_list_of_points ( search         ,
				  true_barrier   ,
				  sgte_barrier   ,
				  NULL           , // Pareto front = NULL
				  stop           ,
				  stop_reason    ,
				  new_feas_inc   ,
				  new_infeas_inc ,
				  success        ,
				  *evaluated_pts   );
    if ( stop ) {
      if ( del_evaluated_pts )
	delete evaluated_pts;
      return;
    }

    NOMAD::Eval_Point * x;

    // construct a new list of trial points that will be
    // ordered using surrogate values:   
    std::list<const NOMAD::Eval_Point *>::const_iterator
      end = evaluated_pts->end() , it2;
    for ( it2 = evaluated_pts->begin() ; it2 != end ; ++it2 ) {
      
      // Eval_Point construction:
      x = new NOMAD::Eval_Point;
      x->set ( (*it2)->size() , _p.get_bb_nb_outputs() );
      x->set_signature  ( (*it2)->get_signature () );
      x->set_direction  ( (*it2)->get_direction () );
      x->set_mesh_index ( (*it2)->get_mesh_index() );
      x->Point::operator = ( **it2 );
	
      display_new_list = true;
	
      // add the new point to the ordered list of trial points:
      add_eval_point ( x                       ,
		       display_degree          ,
		       _p.get_snap_to_bounds() ,
		       (*it2)->get_f()         ,
		       (*it2)->get_h()           );
    }
  }
  
  if ( stop ) {
    if ( del_evaluated_pts )
      delete evaluated_pts;
    return;
  }

  // display the re-ordered list of trial points:
  if ( display_new_list && display_degree == NOMAD::FULL_DISPLAY ) {

    const NOMAD::Eval_Point * y;

    std::ostringstream oss;
    oss << "re-ordered list of " << _eval_lop.size()
	<< " " << search << " trial points";

    out << NOMAD::open_block ( oss.str() ) << std::endl;

    std::set<NOMAD::Priority_Eval_Point>::const_iterator
      end = _eval_lop.end() , it;
    for ( it = _eval_lop.begin() ; it != end ; ++it ) {
      y =  it->get_point();
      y->display_tag ( out );
      out << " : ( ";
      y->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
      out << " )";
      if ( y->get_direction() )
	out << " (dir " << y->get_direction()->get_index() << ")";
      out << std::endl;
    }
    out.close_block();
  }

  // evaluate the list of points on the 'true' function:
  private_eval_list_of_points ( search         ,
				true_barrier   ,
				sgte_barrier   ,
				pareto_front   ,
				stop           ,
				stop_reason    ,
				new_feas_inc   ,
				new_infeas_inc ,
				success        ,
				*evaluated_pts   ); 

  if ( del_evaluated_pts )
    delete evaluated_pts;
}

/*--------------------------------------------------------------*/
/*  return if a series of evaluations is opportunistic or not,  */
/*  depending on the search type (private)                      */
/*--------------------------------------------------------------*/
bool NOMAD::Evaluator_Control::is_opportunistic ( NOMAD::search_type t ) const
{
  switch ( t ) {
  case NOMAD::X0_EVAL:
    return false;
  case NOMAD::LH_SEARCH:
    return _p.get_opportunistic_LH();
  case NOMAD::CACHE_SEARCH:
    return _p.get_opportunistic_cache_search();
  default:
    return _p.get_opportunistic_eval();
  }
  return false;
}

/*----------------------------------------------------------------*/
/*                     stop the evaluations ?                     */
/*----------------------------------------------------------------*/
/* . check the opportunistic strategy stopping criterion          */
/* . private method                                               */
/*----------------------------------------------------------------*/
bool NOMAD::Evaluator_Control::stop_evaluations
( const NOMAD::Eval_Point & x                ,
  NOMAD::search_type        search           ,
  int                       k                ,
  int                       nb_points        ,
  bool                      stop             ,
  NOMAD::dd_type            display_degree   ,
  NOMAD::success_type       one_eval_success ,
  NOMAD::success_type       success          ,
  int                       init_nb_eval     ,
  const NOMAD::Double     & f0               ,
  const NOMAD::Barrier    & barrier          ,
  int                     & nb_success       ,
  bool                    & one_for_luck       ) const
{
  // opportunistic evaluation ?
  bool opportunistic = is_opportunistic ( search );

  if ( k < nb_points - 1 ) {
	
    if ( stop )
      return true;

    if ( opportunistic &&
	 ( x.get_eval_type() == NOMAD::TRUTH || _p.get_opt_only_sgte() ) ) {
      
      if ( one_for_luck && one_eval_success != NOMAD::FULL_SUCCESS ) {
	if ( display_degree == NOMAD::FULL_DISPLAY )
	  _p.out() << std::endl
		   << "opportunistic termination of evaluations (lucky eval)"
		   << std::endl;
	return true;
      }
      
      if ( success == NOMAD::FULL_SUCCESS &&
	   check_opportunistic_criterion ( display_degree   ,
					   one_eval_success ,
					   init_nb_eval     ,
					   f0               ,
					   barrier          ,
					   nb_success       ,
					   one_for_luck       ) )
	return true;
    }
  }
  return false;
}

/*-----------------------------------------------------------------*/
/*  check the opportunistic strategy stopping criterion (private)  */
/*            return true to stop the evaluations                  */
/*            return false to continue the evaluations             */
/*-----------------------------------------------------------------*/
bool NOMAD::Evaluator_Control::check_opportunistic_criterion
( NOMAD::dd_type         display_degree   ,
  NOMAD::success_type    one_eval_success ,
  int                    init_nb_eval     ,
  const NOMAD::Double  & f0               ,
  const NOMAD::Barrier & barrier          ,
  int                  & nb_success       ,
  bool                 & one_for_luck       ) const
{
  int                    min_nb_success = _p.get_opportunistic_min_nb_success();
  int                    min_eval       = _p.get_opportunistic_min_eval();
  NOMAD::Double          min_f_imprvmt  = _p.get_opportunistic_min_f_imprvmt();
  bool                   lucky_eval     = _p.get_opportunistic_lucky_eval();
  const NOMAD::Display & out            = _p.out();

  // min_nb_success:
  if ( min_nb_success > 0 ) {

    if ( one_eval_success == NOMAD::FULL_SUCCESS )
      ++nb_success;

    if ( nb_success < min_nb_success ) {

      if ( display_degree == NOMAD::FULL_DISPLAY )
	out << std::endl
	    << "opport. strategy (nb_success=" << nb_success
	    << " < min_nb_success=" << min_nb_success
	    << "): continue evaluations"
	    << std::endl;
      
      return false;
    }
  }

  // min_eval:
  if ( min_eval > 0 ) {
    
    int eval = _stats.get_eval() - init_nb_eval;

    if ( eval < min_eval ) {

      if ( display_degree == NOMAD::FULL_DISPLAY )
	out << std::endl
	    << "opport. strategy (eval=" << eval
	    << " < min_eval=" << min_eval
	    << "): continue evaluations" << std::endl;
      return false;
    }
  }

  // min_f_imprvmt:
  if ( min_f_imprvmt.is_defined() ) {

    const NOMAD::Eval_Point * bf = barrier.get_best_feasible();

    if ( f0.is_defined() && bf ) {

      NOMAD::Double f = bf->get_f();

      if ( f.is_defined() ) {
		  
	NOMAD::Double f_imprvmt = ( f != 0.0 ) ?
	  100.0 * (f0 - f) / f.abs() : 100.0;
	
	if ( f_imprvmt < min_f_imprvmt ) {

	  if ( display_degree == NOMAD::FULL_DISPLAY )
	    out << std::endl
		<< "opport. strategy (f_improvement="
		<< f_imprvmt << " < min_f_imprvmt=" << min_f_imprvmt
		<< "): continue evaluations" << std::endl;

	  return false;
	}
      }
    }
  }

  // lucky_eval:
  if ( lucky_eval && one_eval_success == NOMAD::FULL_SUCCESS ) {
    one_for_luck = true;

    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << std::endl
	  << "opport. strategy: one more evaluation for luck"
	  << std::endl;

    return false;
  }

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl << "opport. strategy: stop evaluations"
	<< std::endl;

  return true;
}

/*---------------------------------------------------------------*/
/*        display the list of evaluation points (_eval_lop)      */
/*---------------------------------------------------------------*/
void NOMAD::Evaluator_Control::display_eval_lop ( NOMAD::search_type t ) const
{
  const NOMAD::Display & out = _p.out();
  int cnt = 0 , nb = _eval_lop.size();

  if ( nb == 0 ) {
    out << std::endl << "no evaluation point" << std::endl;
    return;
  }

  // open indented block:
  std::ostringstream oss;
  if ( t != NOMAD::UNDEFINED_SEARCH )
    oss << t << " ";
  oss << "evaluation point";
  if ( nb > 1 )
    oss << "s";
  out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;

  // display the points:
  std::set<NOMAD::Priority_Eval_Point>::const_iterator it , end = _eval_lop.end();
  for ( it = _eval_lop.begin() ; it != end ; ++it ) {
    out << "point ";
    out.display_int_w ( ++cnt , nb );
    out << "/" << nb << ": ( ";
    it->get_point()->Point::display ( out                               ,
				      " "                               ,
				      2                                 ,
				      NOMAD::Point::get_display_limit()   );
    out << " )" << std::endl;
  }

  // close indented block:
  out.close_block();
}

/*--------------------------------------------------------------*/
/*    add an Eval_Point to the list of points to be evaluated   */
/*--------------------------------------------------------------*/
/*  . x has to be a dynamic object                              */
/*  . it may be deleted into the method and be NULL after that  */
/*  . the point is also snapped to bounds                       */
/*  . periodic variables are checked                            */
/*--------------------------------------------------------------*/
void NOMAD::Evaluator_Control::add_eval_point
( NOMAD::Eval_Point  *& x              ,
  NOMAD::dd_type        display_degree ,
  bool                  snap_to_bounds ,
  const NOMAD::Double & f_sgte         ,
  const NOMAD::Double & h_sgte           )
{
  if ( !x )
    return;

  const NOMAD::Display & out = _p.out();

  // treat the periodic variables:
  NOMAD::Direction * new_dir = NULL;
  
  if ( _p.has_periodic_variables() &&
       x->treat_periodic_variables ( new_dir ) ) {

    if ( new_dir && new_dir->norm() == 0.0 ) {
      
      if ( display_degree == NOMAD::FULL_DISPLAY )
	out << "point #" << x->get_tag()
	    << " is flushed (||dir||==0)"
	    << std::endl;

      delete x;
      x = NULL;

      delete new_dir;

      return;
    }
  }
  delete new_dir;

  // snap to bounds:
  if ( snap_to_bounds && x->snap_to_bounds() ) {

    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      out << std::endl << "point #" << x->get_tag() << " ";
      if ( x->get_direction() && x->get_direction()->get_index() >= 0 )
	out << "(dir " << x->get_direction()->get_index() << ") ";
      out << "has been snapped to bounds" << std::endl;
    }

    if ( x->get_direction() && x->get_direction()->norm() == 0.0 ) {
      
      if (display_degree == NOMAD::FULL_DISPLAY )
	out << "point #" << x->get_tag()
	    << " is flushed (||dir||==0)"
	    << std::endl;
      delete x;
      x = NULL;

      return;
    }
  }

  // creation of the Priority_Eval_Point:
  NOMAD::Priority_Eval_Point pep ( x , _p.get_h_min() );

  // ordering elements of Priority_Eval_Point's:
  // -------------------------------------------

  // 1. surrogate values for f and h:
  pep.set_f_sgte ( f_sgte );
  pep.set_h_sgte ( h_sgte );
  
  // 2. angle with last successful directions:
  if ( x->get_direction() ) {
    
    // get the signature:
    NOMAD::Signature * signature = x->get_signature();
    if ( !signature )
      throw NOMAD::Exception ( "Evaluator_Control.cpp" , __LINE__ ,
	    "Evaluator_Control::add_eval_point(): the point has no signature" );
    
    // feasible success direction:
    const NOMAD::Direction & feas_success_dir = signature->get_feas_success_dir();
    if ( feas_success_dir.is_defined() &&
	 x->get_poll_center_type() == NOMAD::FEASIBLE )
      pep.set_angle_success_dir ( feas_success_dir.get_angle( *x->get_direction() ) );

    // infeasible success direction:
    const NOMAD::Direction & infeas_success_dir = signature->get_infeas_success_dir();
    if ( infeas_success_dir.is_defined() &&
	 x->get_poll_center_type() == NOMAD::INFEASIBLE )
      pep.set_angle_success_dir ( infeas_success_dir.get_angle( *x->get_direction() ) );
  }

  // 3. angle with simplex gradient:
  // (to be implemented)

  // insertion of the point in _eval_lop:
  // ------------------------------------
  _eval_lop.insert ( pep );
}
