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
  \file   LH_Search.cpp
  \brief  Latin-Hypercube search (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    LH_Search.hpp
*/
#include "LH_Search.hpp"

/*-----------------------------------------------------------*/
/*              MADS Latin-Hypercube (LH) search             */
/*-----------------------------------------------------------*/
void NOMAD::LH_Search::search ( NOMAD::Mads              & mads           ,
				int                      & nb_search_pts  ,
				bool                     & stop           ,
				NOMAD::stop_type         & stop_reason    ,
				NOMAD::success_type      & success        ,
				bool                     & count_search   ,
				const NOMAD::Eval_Point *& new_feas_inc   ,
				const NOMAD::Eval_Point *& new_infeas_inc   )
{
  new_feas_inc = new_infeas_inc = NULL;
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
    oss << NOMAD::LH_SEARCH;
    out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
  }

  // active barrier:
  const NOMAD::Barrier     & barrier = mads.get_active_barrier();

  // Evaluator_Control:
  NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();

  // current incumbents:
  const NOMAD::Eval_Point  * feas_inc   = barrier.get_best_feasible  ();
  const NOMAD::Eval_Point  * infeas_inc = barrier.get_best_infeasible();

  // get a reference point and a signature:
  const NOMAD::Eval_Point  * ref       = (feas_inc) ? feas_inc : infeas_inc;
  NOMAD::Signature         * signature = _p.get_signature();

  // check the number of points:
  int p = _initial_search ? _p.get_LH_search_p0() : _p.get_LH_search_pi();
  if ( p <= 0 ) {
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      std::ostringstream oss;
      oss << "end of LH " << ( _initial_search ? "initial " : "")
	  << "search (number of points <= 0)";
      out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
    }
    return;
  }

  // no reference point is available (we consider the standard signature):
  if ( !ref ) {

    // it is not sufficient with categorical variables:
    if ( signature->has_categorical() ) {
      if ( display_degree == NOMAD::FULL_DISPLAY ) {
	std::ostringstream oss;
	oss << "end of LH " << ( _initial_search ? "initial " : "")
	    << "search (no available reference point)";
	out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
      }
      return;
    }
  }
  else
    signature = ref->get_signature();

  int                      i;
  NOMAD::Eval_Point      * x;
  int                      n          = signature->get_n();
  int                      m          = _p.get_bb_nb_outputs();
  int                      mesh_index = NOMAD::Mesh::get_mesh_index();
  int                      pm1        = p-1;
  
  // mesh sizes:
  NOMAD::Point delta_m_max ( n );
  signature->get_mesh().get_delta_m ( delta_m_max , NOMAD::Mesh::get_min_mesh_index() );
  
  NOMAD::Double delta_m_i;
  NOMAD::Point  delta_m;
  if ( !_initial_search )
    signature->get_mesh().get_delta_m ( delta_m , mesh_index );

  // fixed variables:
  const NOMAD::Point & fixed_variables = signature->get_fixed_variables();
  
  // bb input types:
  const std::vector<NOMAD::bb_input_type> & bbit = signature->get_input_types();
  
  // bounds:
  const NOMAD::Point & lb = signature->get_lb();
  const NOMAD::Point & ub = signature->get_ub();

  // pts contains n points of dimension p: each of these points contains
  // p different values for each variable:
  NOMAD::Point ** pts = new NOMAD::Point * [n];

  // creation of p search points:
  for ( int k = 0 ; k < p ; ++k ) {
    
    x = new NOMAD::Eval_Point ( n , m );
    x->set_signature  ( signature   );
    x->set_mesh_index ( &mesh_index );
    
    for ( i = 0 ; i < n ; ++i ) {
      
      if ( k==0 ) {
	if ( fixed_variables[i].is_defined() )
	  pts[i] = new NOMAD::Point ( p , fixed_variables[i] );
	else if ( bbit[i] == NOMAD::CATEGORICAL ) {
	  pts[i] = new NOMAD::Point ( p , (*ref)[i] );
	}
	else {
	  pts[i] = new NOMAD::Point ( p );
	  
	  // for the initial mesh: delta_m is not used and there will
	  // be no projection on mesh:
	  if ( !_initial_search )
	    delta_m_i = delta_m[i];

	  values_for_var_i ( p               ,
			     delta_m_i       ,
			     delta_m_max[i]  ,
			     bbit       [i]  ,
			     lb         [i]  ,
			     ub         [i]  ,
			     *pts       [i]    );
	}
      }
      
      (*x)[i] = (*pts[i])[k];
      
      if ( k == pm1 )
	delete pts[i];
    }
    
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      out << "LH point #" << x->get_tag()
	  << ": ( ";
      x->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
      out << " )" << std::endl;
    }

    // add the new point to the ordered list of search trial points:
    ev_control.add_eval_point ( x               ,
				display_degree  ,
				false           ,  // snap_to_bounds = false
				NOMAD::Double() ,
				NOMAD::Double()   );
  }
  
  delete [] pts;
  
  nb_search_pts = static_cast<int> ( ev_control.get_eval_lop().size() );

  // eval_list_of_points:
  // --------------------
  new_feas_inc = new_infeas_inc = NULL;
  ev_control.eval_list_of_points ( _type                   ,
				   mads.get_true_barrier() ,
				   mads.get_sgte_barrier() ,
				   mads.get_pareto_front() ,
				   stop                    ,
				   stop_reason             ,
				   new_feas_inc            ,
				   new_infeas_inc          ,
				   success                   );

  // final displays:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {
    std::ostringstream oss;
    oss << "end of LH search (" << success << ")";
    out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
  }

}

/*-----------------------------------------------------------*/
/*        LH search: decide p values for one variable        */
/*-----------------------------------------------------------*/
/*  . if no bounds, values are scaled with the largest       */
/*    delta_m value obtained so far                          */
/*  . private method                                         */
/*-----------------------------------------------------------*/
void NOMAD::LH_Search::values_for_var_i ( int                          p           ,
					  const NOMAD::Double        & delta_m     ,
					  const NOMAD::Double        & delta_m_max ,
					  const NOMAD::bb_input_type & bbit        ,
					  const NOMAD::Double        & lb          ,
					  const NOMAD::Double        & ub          ,
					  NOMAD::Point               & x      ) const
{
  // categorical variables have already been treated as fixed variables:
  if ( bbit == NOMAD::CATEGORICAL )
    return;

  int                  i;
  NOMAD::Double        v;
  NOMAD::Random_Pickup rp ( p );
  bool                 rounding = ( bbit != NOMAD::CONTINUOUS );
  bool                 lb_def   = lb.is_defined();
  bool                 ub_def   = ub.is_defined();
  double               w        = ( ( lb_def && ub_def ) ?
				    ub.value()-lb.value() : 1.0 ) / p;
  // main loop:
  for ( i = 0 ; i < p ; ++i ) {

    // both bounds exist:
    if ( lb_def && ub_def )
      v = lb + ( i + rand()/NOMAD::D_INT_MAX ) * w;

    // one of the bounds does not exist:
    else {

      // lb exists, and ub not: mapping [0;1] --> [lb;+INF[
      if ( lb_def )
	v = lb + 10 *
	  delta_m_max * sqrt ( - log ( NOMAD::DEFAULT_EPSILON +
				       ( i + rand()/NOMAD::D_INT_MAX ) * w ) );

      // lb does not exist:
      else {

	// ub exists, and lb not: mapping [0;1] --> ]-INF;ub]
	if ( ub_def )
	  v = ub - delta_m_max * 10 *
	    sqrt ( -log ( NOMAD::DEFAULT_EPSILON +
			  ( i + rand()/NOMAD::D_INT_MAX ) * w ) );
	
	// there are no bounds: mapping [0;1] --> ]-INF;+INF[
	else
	  v = (rand()%2 ? -1.0 : 1.0) * delta_m_max * 10 *
	    sqrt ( - log ( NOMAD::DEFAULT_EPSILON +
			   ( i + rand()/NOMAD::D_INT_MAX ) * w ) );
      }
    }

    // rounding:
    if ( rounding )
      v = v.round();
   
    // projection to mesh (with ref=0):
    v.project_to_mesh ( 0.0 , delta_m , lb , ub );

    // affectation + permutation:
    x[rp.pickup()] = v;
  }
}
