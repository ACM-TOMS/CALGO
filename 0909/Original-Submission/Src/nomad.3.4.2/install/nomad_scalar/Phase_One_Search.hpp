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
  \file   Phase_One_Search.hpp
  \brief  NOMAD::Search subclass for the phase one (headers)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    Phase_One_Search.cpp
*/
#ifndef __PHASE_ONE_SEARCH__
#define __PHASE_ONE_SEARCH__

#include "Search.hpp"
#include "Mads.hpp"

namespace NOMAD {

  /// NOMAD::Search subclass for the phase one.
  /**
     - The phase one occurs when no feasible starting point has been given.
     - It consists in minimizing the constraint violations and it stops
       as soon as a feasible point is found.
  */
  class Phase_One_Search : public NOMAD::Search , private NOMAD::Uncopyable {

  public:

    /// Constructor.
    /**
       \param p Parameters -- \b IN.
    */
    Phase_One_Search ( NOMAD::Parameters & p )
      : NOMAD::Search ( p , NOMAD::P1_SEARCH ) {}
  
    /// Destructor.
    virtual ~Phase_One_Search ( void ) {}

    /// The phase one search.
    /**
       \param mads           NOMAD::Mads object invoking this search -- \b IN/OUT.
       \param nb_search_pts  Number of generated search points       -- \b OUT.
       \param stop           Stop flag                               -- \b IN/OUT.
       \param stop_reason    Stop reason                             -- \b OUT.
       \param success        Type of success                         -- \b OUT.
       \param count_search   Count or not the search                 -- \b OUT.
       \param new_feas_inc   New feasible incumbent                  -- \b IN/OUT.
       \param new_infeas_inc New infeasible incumbent                -- \b IN/OUT.
    */
    virtual void search ( NOMAD::Mads              & mads           ,
			  int                      & nb_search_pts  ,
			  bool                     & stop           ,
			  NOMAD::stop_type         & stop_reason    ,
			  NOMAD::success_type      & success        ,
			  bool                     & count_search   ,
			  const NOMAD::Eval_Point *& new_feas_inc   ,
			  const NOMAD::Eval_Point *& new_infeas_inc   );
  };
}

#endif
