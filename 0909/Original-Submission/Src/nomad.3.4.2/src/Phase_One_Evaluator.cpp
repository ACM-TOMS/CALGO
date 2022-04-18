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
  \file   Phase_One_Evaluator.cpp
  \brief  NOMAD::Evaluator subclass for the phase one (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    Phase_One_Evaluator.hpp
*/
#include "Phase_One_Evaluator.hpp"

/*------------------------------------------------------------------*/
/*         compute f(x) from the blackbox outputs of a point        */
/*              (special objective for MADS phase 1)                */
/*------------------------------------------------------------------*/
void NOMAD::Phase_One_Evaluator::compute_f ( NOMAD::Eval_Point & x ) const
{
  if ( x.get_bb_outputs().size() != _p.get_bb_nb_outputs() ) {
    std::ostringstream err;
    err << "Phase_One_Evaluator::compute_f(x): "
	<< "x has a wrong number of blackbox outputs ("
	<< x.get_bb_outputs().size() <<  " != " << _p.get_bb_nb_outputs() << ")";
    throw NOMAD::Exception ( "Phase_One_Evaluator.cpp" , __LINE__ , err.str() );
  }

  // objective value for MADS phase 1: the squared sum of all EB constraint violations
  // (each EB constraint has been previously transformed into OBJ values):
  const std::list<int>               & index_obj = _p.get_index_obj();
  const std::list<int>::const_iterator end       = index_obj.end();
  const NOMAD::Point                 & bbo       = x.get_bb_outputs();
  NOMAD::Double                        h_min     = _p.get_h_min();
  NOMAD::Double                        sum       = 0.0;
  NOMAD::Double                        v;

  for ( std::list<int>::const_iterator it = index_obj.begin() ; it != end ; ++it ) {
    v = bbo[*it];
    if ( v > h_min )
      sum += v.pow2();
  }

  x.set_f ( sum );
}
