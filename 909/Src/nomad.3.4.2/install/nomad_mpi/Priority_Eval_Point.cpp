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
  \file   Priority_Eval_Point.cpp
  \brief  Evaluation point with a priority (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-22
  \see    Priority_Eval_Point.hpp
*/
#include "Priority_Eval_Point.hpp"

/*----------------------------------------------------------------*/
/*  . specific Priority_Eval_Point elements of comparison         */
/*  . this method is defined virtual in Set_Element so that       */
/*    Priority_Eval_Point::operator < (Set_Element x) can invoke  */
/*    it on x (which is in fact a Priority_Eval_Point)            */
/*  . this avoids an expensive downcast in operator <             */
/*----------------------------------------------------------------*/
void NOMAD::Priority_Eval_Point::get_priority_criteria
( NOMAD::Double & f_sgte             ,
  NOMAD::Double & h_sgte             ,
  NOMAD::Double & angle_success_dir  ,
  NOMAD::Double & angle_simplex_grad   ) const
{
  f_sgte             = _f_sgte;
  h_sgte             = _h_sgte;
  angle_success_dir  = _angle_success_dir;
  angle_simplex_grad = _angle_simplex_grad;
}

/*------------------------------------------------------*/
/*                  comparison operator                 */
/*  . x1 < x2 true if x1 should be evaluated before x2  */
/*  . x is a Priority_Eval_Point                        */
/*------------------------------------------------------*/
bool NOMAD::Priority_Eval_Point::operator <
( const NOMAD::Set_Element<NOMAD::Eval_Point> & x ) const
{
  if ( this == &x )
    return false;

  const NOMAD::Eval_Point * x1 = get_element();
  const NOMAD::Eval_Point * x2 = x.get_element();
  
  // criterion 1: user criterion:
  // ------------
  const NOMAD::Double uep1 = x1->get_user_eval_priority();
  if ( uep1.is_defined() ) {
    const NOMAD::Double uep2 = x2->get_user_eval_priority();
    if ( uep2.is_defined() ) {
      if ( uep1 > uep2 )
	return true;
      if ( uep2 > uep1 )
	return false;
    }
  }

  // specific Priority_Eval_Point elements of comparison:
  NOMAD::Double x_f_sgte;             
  NOMAD::Double x_h_sgte;             
  NOMAD::Double x_angle_success_dir;
  NOMAD::Double x_angle_simplex_grad;

  x.get_priority_criteria ( x_f_sgte             ,
			    x_h_sgte             ,
			    x_angle_success_dir  ,
			    x_angle_simplex_grad   );

  // criterion 2: give priority to already evaluated cache points:
  // ------------
  if ( x1->is_in_cache() && !x2->is_in_cache() )
    return true;
  if ( x2->is_in_cache() && !x1->is_in_cache() )
    return false;

  // criterion 3: give priority to already evaluated points
  // ------------ that are eval_ok:
  if ( x1->is_eval_ok() && !x2->is_eval_ok() )
    return true;
  if ( x2->is_eval_ok() && !x1->is_eval_ok() )
    return false;

  // criterion 4: true f and h values:
  // -----------
  int flag = compare_hf_values ( x1->get_h() ,
				 x1->get_f() ,
				 x2->get_h() ,
				 x2->get_f()   );
  if ( flag )
    return ( flag > 0 );

  // criterion 5: surrogate f and h values:
  // ------------
  flag = compare_hf_values ( _h_sgte , _f_sgte , x_h_sgte , x_f_sgte );
  if ( flag )
    return ( flag > 0 );

  // criterion 6: check the angle with the last successful direction:
  // ------------
  if ( _angle_success_dir.is_defined() && x_angle_success_dir.is_defined() ) {

    if ( _angle_success_dir < x_angle_success_dir )
      return true;
  
    if ( x_angle_success_dir < _angle_success_dir )
      return false;
  }

  // criterion 7: check the angle with the simplex gradient:
  // ------------
  if ( _angle_simplex_grad.is_defined() && x_angle_simplex_grad.is_defined() ) {
    
    if ( _angle_simplex_grad < x_angle_simplex_grad )
      return true;

    if ( x_angle_simplex_grad < _angle_simplex_grad )
      return false;
  }

  // criterion 8: random criterion for randomly generated directions:
  // ------------
  const NOMAD::Double rep1 = x1->get_rand_eval_priority();
  if ( rep1.is_defined() ) {
    const NOMAD::Double rep2 = x2->get_rand_eval_priority();
    if ( rep2.is_defined() ) {
      if ( rep1 < rep2 )
	return true;
      if ( rep2 < rep1 )
	return false;
    }
  }

  // criterion 9: compare the tags:
  // ------------
  return x1->get_tag() < x2->get_tag();
}

/*-----------------------------------------------*/
/*    compare the h and f values of two points   */
/*-----------------------------------------------*/
/*  . return ( h(x1),f(x1) ) < ( h(x2),f(x2) )   */
/*    with the following format:                 */
/*                      1 : x1 best than x2      */
/*                     -1 : x2 best than x1      */
/*                      0 : undetermined         */
/*  . private method                             */
/*-----------------------------------------------*/
int NOMAD::Priority_Eval_Point::compare_hf_values ( const NOMAD::Double & hx1 ,
						    const NOMAD::Double & fx1 ,
						    const NOMAD::Double & hx2 ,
						    const NOMAD::Double & fx2   ) const
{
  if ( fx1.is_defined() && fx2.is_defined() ) {

    if ( hx1.is_defined() && hx2.is_defined() ) {

      // x1 is feasible:
      if ( hx1 <= _h_min ) {

	// both points are feasible:
	if ( hx2 <= _h_min  ) {
	  if ( fx1 < fx2 )
	    return 1;
	  if ( fx2 < fx1 )
	    return -1;
	}
	
	// x1 feasible and x2 infeasible:
	else
	  return 1;
      }

      // x1 is infeasible:
      else {
	
	// x2 is feasible:
	if ( hx2 <= _h_min  )
	  return -1;
	
	// both points are infeasible:
	if ( ( hx1  < hx2 && fx1  < fx2 ) ||
	     ( hx1 == hx2 && fx1  < fx2 ) ||
	     ( hx1  < hx2 && fx1 == fx2 )    )
	  return 1;

	if ( ( hx2  < hx1 && fx2  < fx1 ) ||
	     ( hx2 == hx1 && fx2  < fx1 ) ||
	     ( hx2  < hx1 && fx2 == fx1 )    )
	  return -1; 
      }
    }

    // we only have f values:
    else {
      if ( fx1 < fx2 )
	return 1;
      if ( fx2 < fx1 )
	return -1;
    }
  }

  return 0;
}
