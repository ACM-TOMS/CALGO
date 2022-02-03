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
  \file   VNS_Search.cpp
  \brief  VNS search (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-12
  \see    VNS_Search.hpp
*/
#include "VNS_Search.hpp"

/*---------------------------------------------------------*/
/*                           reset                         */
/*---------------------------------------------------------*/
void NOMAD::VNS_Search::reset ( void )
{
  _k = _k_max   = 1;
  _halton_index = NOMAD::VNS_HALTON_INDEX_0;
  _old_x        = NULL;
  for ( int ell = 0 ; ell <= NOMAD::L_LIMITS ; ++ell )
    _nb_performed[ell] = 0;
}

/*---------------------------------------------------------*/
/*                       the search                        */
/*       VNS: x --[shaking(k)]--> x' --[descent]--> x"     */
/*---------------------------------------------------------*/
void NOMAD::VNS_Search::search ( NOMAD::Mads              & mads           ,
				 int                      & nb_search_pts  ,
				 bool                     & stop           ,
				 NOMAD::stop_type         & stop_reason    ,
				 NOMAD::success_type      & success        ,
				 bool                     & count_search   ,
				 const NOMAD::Eval_Point *& new_feas_inc   ,
				 const NOMAD::Eval_Point *& new_infeas_inc   )
{
  new_feas_inc  = new_infeas_inc = NULL;
  nb_search_pts = 0;
  success       = NOMAD::UNSUCCESSFUL;
  count_search  = !stop;

  if ( stop )
    return;

  // initial display:
  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_search_dd();
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    std::ostringstream oss;
    oss << NOMAD::VNS_SEARCH;
    out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
  }

  bool opt_only_sgte = _p.get_opt_only_sgte();

  // the barriers:
  NOMAD::Barrier       & true_barrier   = mads.get_true_barrier();
  NOMAD::Barrier       & sgte_barrier   = mads.get_sgte_barrier();
  const NOMAD::Barrier & active_barrier = mads.get_active_barrier();

  // point x:
  NOMAD::Double             best_f;
  bool                      x_feas = true;
  const NOMAD::Eval_Point * x      = active_barrier.get_best_feasible();
  if ( x )
    best_f = x->get_f();
  else {
    x      = active_barrier.get_best_infeasible();
    x_feas = false;
  }

  if ( !x ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out.close_block ( "end of VNS search (no incumbent)" );
    return;
  }

  // update _k and _old_x:
  if ( x == _old_x ) {
    ++_k;
    if ( _k > _k_max )
      _k_max = _k;
  }
  else
    _k = 1;
  
  _old_x = x;

  // get the signature:
  NOMAD::Signature * signature = x->get_signature();
  if ( !signature )
    return;

  int n = signature->get_n();
  if ( n != x->size() ) {
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out.close_block ( "end of VNS search (incompatible signature)" );
    return;
  }

  // shaking: get ONE direction from the signature:
  NOMAD::Direction dir;
  signature->get_one_direction ( dir           ,
				 1 - _k        ,   // mesh_index = 1-k
				 _halton_index   );

  _halton_index += NOMAD::VNS_HALTON_INCR;
  
  // shaking: construct x':
  NOMAD::Point xp = *x + dir;

  // shaking: the perturbation is tried twice with dir and -dir
  //          (in case x == x + dir after snapping)
  for ( int nbt = 0 ; nbt < 2 ; ++nbt ) {

    // treat xp: periodic variables or bounds:
    if ( _p.has_periodic_variables() ) {
      NOMAD::Direction * tmp_dir = NULL;
      signature->treat_periodic_variables ( xp , NULL , tmp_dir );
    }
    else
      signature->snap_to_bounds ( xp , NULL );

    // if xp == x :
    if ( xp == *x ) {
      
      // no third try: the search fails
      if ( nbt == 1 ) {
	if ( display_degree == NOMAD::FULL_DISPLAY )
	  out.close_block ( "end of VNS search (shaking failed)" );
	return;
      }
      
      // 2nd try (-dir instead of dir):
      xp = *x - dir;
    }
  }

  // current mesh index:
  int mesh_index = NOMAD::Mesh::get_mesh_index();

  // reset mesh index:
  int initial_mesh_index = _p.get_initial_mesh_index();
  NOMAD::Mesh::set_mesh_index ( initial_mesh_index );

  // stats:
  NOMAD::Stats & stats = mads.get_stats();

  // current number of blackbox evaluations:
  int  bbe             = stats.get_bb_eval();
  int  sgte_eval       = stats.get_sgte_eval();
  int  mads_iterations = stats.get_iterations();
  bool has_sgte        = _p.has_sgte();

  // displays:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    out << "        it = " << mads_iterations << std::endl
	<< "       bbe = " << bbe << std::endl;
    if ( has_sgte )
      out << " sgte_eval = " << sgte_eval << std::endl;
    out << "mesh_index = " << mesh_index << std::endl
	<< "         k = " << _k << std::endl
	<< "      kmax = " << _k_max << std::endl
	<< "         x = ( ";
    x->Point::display ( out , " " , 5 , _p.get_point_display_limit() );
    out << " ) f=" << x->get_f() << " h=" << x->get_h() << std::endl
	<< "       dir = ( ";
    dir.Point::display ( out , " " , 5 , _p.get_point_display_limit() );
    out << " ) |dir|=";
    NOMAD::Double norm = dir.norm();
    out << norm << std::endl; 
    out << "        xp = ( ";
    xp.display ( out , " " , 5 , _p.get_point_display_limit() );
    out << " )" << std::endl << std::endl;
    out << "bb_eval (before+VNS only) objective_value"
	<< std::endl << std::endl;
  }

  // save parameters that are going to be modified:
  // ----------------------------------------------
  std::string old_display_degree;
  _p.out().get_display_degree ( old_display_degree );
  bool                                old_ses = _p.get_sgte_eval_sort();
  bool                                old_sif = _p.get_stop_if_feasible();
  int                                  old_hs = _p.get_halton_seed();
  int                            old_max_time = _p.get_max_time();
  int                      old_max_mesh_index = _p.get_max_mesh_index();
  int                             old_max_bbe = _p.get_max_bb_eval();
  int                            old_max_eval = _p.get_max_eval();
  int                       old_max_sgte_eval = _p.get_max_sgte_eval();
  int                              old_max_it = _p.get_max_iterations();
  int                               old_LH_p0 = _p.get_LH_search_p0();
  int                               old_LH_pi = _p.get_LH_search_pi();
  bool                             old_opp_LH = _p.get_opportunistic_LH();
  bool                                 old_CS = _p.get_cache_search();
  bool                             old_opp_CS = _p.get_opportunistic_cache_search();
  int                             old_max_sbe = _p.get_max_sim_bb_eval();
  NOMAD::Double                       old_sst = _p.get_stat_sum_target();
  NOMAD::Double                       old_lct = _p.get_L_curve_target();
  NOMAD::Double                   old_trigger = _p.get_VNS_trigger();
  NOMAD::Point                         old_ft = _p.get_f_target();
  const std::list<std::string>         old_ds = _p.get_display_stats();
  const std::list<std::string> old_stats_file = _p.get_stats_file();
  const std::string       old_stats_file_name = _p.get_stats_file_name();
  const std::string              old_sol_file = _p.get_solution_file();
  const std::string              old_his_file = _p.get_history_file();
  bool                                old_uce = _p.get_user_calls_enabled();
  bool                                old_epe = _p.get_extended_poll_enabled();

  // save list of starting points:
  std::string x0_cache_file = _p.get_x0_cache_file();
  std::vector<NOMAD::Point *> x0s;
  {
    const std::vector<NOMAD::Point *> & x0s_tmp = _p.get_x0s();
    size_t nx0 = x0s_tmp.size() , k;
    for ( k = 0 ; k < nx0 ; ++k )
      x0s.push_back ( new Point ( *x0s_tmp[k] ) );
  }

  // modify parameters:
  // ------------------
  _p.set_DISPLAY_DEGREE ( ( display_degree == NOMAD::FULL_DISPLAY ) ?
			  NOMAD::NORMAL_DISPLAY : NOMAD::NO_DISPLAY );
  _p.set_HALTON_SEED    ( NOMAD::Mesh::get_max_halton_index()       );
  _p.set_SOLUTION_FILE  ( ""       );
  _p.set_HISTORY_FILE   ( ""       );
  _p.set_LH_SEARCH      ( 0 , 0    );
  _p.set_VNS_SEARCH     ( false    );
  _p.set_CACHE_SEARCH   ( false    );
  _p.set_MAX_ITERATIONS ( -1       );
  _p.set_OPT_ONLY_SGTE  ( has_sgte );
  _p.reset_X0();
  _p.reset_stats_file();

  _p.set_USER_CALLS_ENABLED    ( false );
  _p.set_EXTENDED_POLL_ENABLED ( false );

  // DISPLAY_STATS:
  {
    if ( has_sgte )
      _p.set_DISPLAY_STATS ( "SGTE OBJ (VNS--surrogate)" );
    else {
      std::list<std::string>                 ds    = old_ds;
      std::list<std::string>::iterator       it    = ds.begin();
      std::list<std::string>::const_iterator end   = ds.end();
      std::string                            s_bbe = NOMAD::itos(bbe) + "+";
      while ( it != end ) {
	if ( *it == "BBE" )
	  ds.insert ( it , s_bbe );
	++it;
      }
      ds.push_back ( " (VNS)" );
      _p.set_DISPLAY_STATS ( ds );
    }
  }

  // mesh:
  _p.set_MAX_MESH_INDEX ( (mesh_index > initial_mesh_index) ?
			  mesh_index : initial_mesh_index);
 
  // X0:
  _p.set_EXTERN_SIGNATURE ( signature );
  _p.set_X0 ( xp );

  // MAX_BB_EVAL:
  if ( old_max_bbe < 0 )
    _p.set_MAX_BB_EVAL ( 100 * n );
  else
    _p.set_MAX_BB_EVAL ( old_max_bbe - bbe );

  // MAX_SGTE_EVAL:
  if ( old_max_sgte_eval > 0 )
    _p.set_MAX_SGTE_EVAL ( old_max_sgte_eval - sgte_eval );
  
  // MAX_EVAL:
  if ( old_max_eval > 0 )
    _p.set_MAX_EVAL ( old_max_eval - stats.get_eval() );

  // MAX_SIM_BB_EVAL:
  if ( old_max_sbe > 0 )
    _p.set_MAX_SIM_BB_EVAL ( old_max_sbe - stats.get_sim_bb_eval() );

  // STAT_SUM_TARGET:
  if ( old_sst.is_defined() )
    _p.set_STAT_SUM_TARGET ( old_sst - stats.get_stat_sum() );

  // MAX_TIME:
  if ( old_max_time > 0 )
    _p.set_MAX_TIME ( old_max_time - stats.get_real_time() );

  // L_CURVE_TARGET:
  if ( !has_sgte )
    _p.set_L_CURVE_TARGET ( best_f );
  
  // F_TARGET and STOP_IF_FEASIBLE:
  if ( has_sgte ) {
    _p.reset_f_target();
    _p.set_STOP_IF_FEASIBLE ( false );
  }

  // check the parameters:
  _p.check ( false ,    // remove_history_file  = false
	     false ,    // remove_solution_file = false
	     false   ); // remove_stats_file    = false

  // Evaluator_Control:
  NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();

  // descent: run MADS:
  // ------------------
  NOMAD::Mads VNS_mads ( _p                           ,
			 ev_control.get_evaluator  () ,
			 mads.get_extended_poll    () ,
			 &ev_control.get_cache     () ,
			 &ev_control.get_sgte_cache()   );

  NOMAD::Mads::set_flag_reset_mesh     ( false );
  NOMAD::Mads::set_flag_reset_barriers ( true  );

  NOMAD::stop_type st = VNS_mads.run();

  NOMAD::Mads::set_flag_reset_mesh ( true );

  // update stats:
  {
    const NOMAD::Stats & VNS_stats = VNS_mads.get_stats();
    stats.update            ( VNS_stats , true ); // for_search = true
    stats.add_VNS_bb_eval   ( VNS_stats.get_bb_eval  () );
    stats.add_VNS_sgte_eval ( VNS_stats.get_sgte_eval() );
  }

  // check MADS stopping criteria:
  if ( st == NOMAD::CTRL_C                   ||
       st == NOMAD::ERROR                    ||
       st == NOMAD::UNKNOWN_STOP_REASON      ||
       st == NOMAD::FEAS_REACHED             ||
       st == NOMAD::MAX_CACHE_MEMORY_REACHED ||
       st == NOMAD::STAT_SUM_TARGET_REACHED  ||
       st == NOMAD::MAX_SGTE_EVAL_REACHED    ||
       st == NOMAD::F_TARGET_REACHED         ||
       st == NOMAD::MAX_SIM_BB_EVAL_REACHED  ||
       st == NOMAD::MAX_TIME_REACHED         ||
       (st == NOMAD::MAX_BB_EVAL_REACHED && old_max_bbe > 0 ) ) {
    stop_reason = st;
    stop        = true;
  }

  // Pareto front:
  NOMAD::Pareto_Front * pareto_front = mads.get_pareto_front();

  // restore starting points:
  {
    _p.reset_X0();
    size_t nx0 = x0s.size();
    if ( nx0 > 0 ) {
      for ( size_t k = 0 ; k < nx0 ; ++k ) {
	_p.set_X0 ( *x0s[k] );
	delete x0s[k];
      }
    }
    else
      _p.set_X0 ( x0_cache_file );
  }

  // restore other saved parameters:
  _p.set_USER_CALLS_ENABLED         ( old_uce                              );
  _p.set_EXTENDED_POLL_ENABLED      ( old_epe                              );
  _p.set_VNS_SEARCH                 ( old_trigger                          );
  _p.set_F_TARGET                   ( old_ft                               );
  _p.set_STOP_IF_FEASIBLE           ( old_sif                              );
  _p.set_L_CURVE_TARGET             ( old_lct                              );
  _p.set_DISPLAY_DEGREE             ( old_display_degree                   );
  _p.set_DISPLAY_STATS              ( old_ds                               );
  _p.set_STATS_FILE                 ( old_stats_file_name , old_stats_file );
  _p.set_SOLUTION_FILE              ( old_sol_file                         );
  _p.set_HISTORY_FILE               ( old_his_file                         );
  _p.set_MAX_BB_EVAL                ( old_max_bbe                          );
  _p.set_MAX_EVAL                   ( old_max_eval                         );
  _p.set_MAX_SGTE_EVAL              ( old_max_sgte_eval                    );
  _p.set_MAX_ITERATIONS             ( old_max_it                           );
  _p.set_STAT_SUM_TARGET            ( old_sst                              );
  _p.set_LH_SEARCH                  ( old_LH_p0 , old_LH_pi                );
  _p.set_OPPORTUNISTIC_LH           ( old_opp_LH                           );
  _p.set_CACHE_SEARCH               ( old_CS                               );
  _p.set_OPPORTUNISTIC_CACHE_SEARCH ( old_opp_CS                           );
  _p.set_MAX_SIM_BB_EVAL            ( old_max_sbe                          );
  _p.set_MAX_MESH_INDEX             ( old_max_mesh_index                   );
  _p.set_HALTON_SEED                ( old_hs                               );
  _p.set_MAX_TIME                   ( old_max_time                         );
  _p.set_SGTE_EVAL_SORT             ( old_ses                              );
  _p.set_OPT_ONLY_SGTE              ( opt_only_sgte                        );
  
  _p.check ( false ,    // remove_history_file  = false
	     false ,    // remove_solution_file = false
	     false   ); // remove_stats_file    = false
    
  // restore mesh index:
  NOMAD::Mesh::set_mesh_index ( mesh_index ); 
  
  // surrogate evaluations: perform only one true evaluation:
  if ( has_sgte && !opt_only_sgte ) {

    if ( !stop ) {

      // remember old best surrogates incumbents:
      const NOMAD::Eval_Point * old_sgte_bf = sgte_barrier.get_best_feasible  ();
      const NOMAD::Eval_Point * old_sgte_bi = sgte_barrier.get_best_infeasible();

      // update the surrogate barrier
      // (no need to invoke Evaluator_Control::process_barrier_points() here
      //  since only surrogate evaluations have been made):
      sgte_barrier.insert ( VNS_mads.get_sgte_barrier() );
      NOMAD::success_type sgte_succ = sgte_barrier.get_success();
      sgte_barrier.update_and_reset_success();

      // we generate only a true trial point if the
      // surrogates improved the surrogate barrier:
      if ( sgte_succ != NOMAD::UNSUCCESSFUL ) {
      
	// choose the best surrogate point(s) where to evaluate the true function:
	const NOMAD::Eval_Point * sgte_bf = sgte_barrier.get_best_feasible  ();
	const NOMAD::Eval_Point * sgte_bi = sgte_barrier.get_best_infeasible();

	std::list<const NOMAD::Eval_Point *> candidates;

	if ( sgte_bf && ( !x_feas || sgte_bf != old_sgte_bf ) )
	  candidates.push_back ( sgte_bf );
      
	if ( sgte_bi && sgte_bi != old_sgte_bi )
	  candidates.push_back ( sgte_bi );

	// generate the new trial points:
	NOMAD::Eval_Point * sk;
	std::list<const NOMAD::Eval_Point *>::const_iterator
	  it , end = candidates.end();
	for ( it = candidates.begin() ; it != end ; ++it ) {

	  // display:
	  if ( display_degree == NOMAD::FULL_DISPLAY )
	    out << std::endl << "VNS surrogate candidate: "
		<< **it << std::endl;

	  sk = new NOMAD::Eval_Point;
	  sk->set ( n , _p.get_bb_nb_outputs() );
	  sk->set_signature  ( signature   );
	  sk->set_mesh_index ( &mesh_index );
      	  sk->Point::operator = ( **it );
	  
	  // add the new point to the list of search trial points:
	  ev_control.add_eval_point ( sk                      ,
				      display_degree          ,
				      _p.get_snap_to_bounds() ,
				      NOMAD::Double()                ,
				      NOMAD::Double()                  );
	}

 	// eval_list_of_points:
 	// --------------------
	success = NOMAD::UNSUCCESSFUL;
 	new_feas_inc = new_infeas_inc = NULL;
	  
 	ev_control.eval_list_of_points ( _type          ,
					 true_barrier   ,
					 sgte_barrier   ,
					 pareto_front   ,
					 stop           ,
					 stop_reason    ,
					 new_feas_inc   ,
					 new_infeas_inc ,
					 success          );

	// number of search points (0 or 1 or 2):
	nb_search_pts = static_cast<int> ( candidates.size() );
      }
    }
  }

  // true evaluations (or surrogate
  // evaluations if opt_only_sgte==true):
  else {

    // for the update of new_feas_inc and new_infeas_inc (1/2):
    const NOMAD::Eval_Point * old_feasible_incumbent   =
      active_barrier.get_best_feasible();
    const NOMAD::Eval_Point * old_infeasible_incumbent =
      active_barrier.get_best_infeasible();
  
    // update barriers and process VNS search points:
    NOMAD::success_type sgte_succ
      = ev_control.process_barrier_points ( sgte_barrier                ,
					    VNS_mads.get_sgte_barrier() ,
					    pareto_front                ,
					    display_degree              ,
					    NOMAD::VNS_SEARCH             );
    NOMAD::success_type true_succ
      = ev_control.process_barrier_points ( true_barrier                ,
					    VNS_mads.get_true_barrier() ,
					    pareto_front                ,
					    display_degree              ,
					    NOMAD::VNS_SEARCH             );
    
    // update of new_feas_inc and new_infeas_inc (2/2):
    const NOMAD::Eval_Point * bf = active_barrier.get_best_feasible  ();
    const NOMAD::Eval_Point * bi = active_barrier.get_best_infeasible();
    if ( bf && bf != old_feasible_incumbent )
      new_feas_inc = bf;
    if ( bi && bi != old_infeasible_incumbent )
      new_infeas_inc = bi;

    // number of search points and success:
    if ( opt_only_sgte ) {
      nb_search_pts = VNS_mads.get_stats().get_sgte_eval();
      success       = sgte_succ;
    }
    else {
      nb_search_pts = VNS_mads.get_stats().get_eval();
      success       = true_succ;
    }
  }

  // update _nb_performed:
  if ( mesh_index > 0 )
    ++_nb_performed [ mesh_index ];

  // final display:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    std::ostringstream oss;
    oss << "end of VNS search (" << success << ")";
    out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
  }
}
