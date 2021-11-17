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
  \file   Cache_Search.hpp
  \brief  NOMAD::Search subclass for the cache search (headers)
  \author Sebastien Le Digabel
  \date   2010-04-08
  \see    Cache_Search.cpp
*/
#ifndef __CACHE_SEARCH__
#define __CACHE_SEARCH__

#include "Mads.hpp"
#include "Search.hpp"

namespace NOMAD {

  /// Cache search.
  class Cache_Search : public NOMAD::Search , private NOMAD::Uncopyable {

  private:

    /// Number of extern points at the end of the last cache search.
    int _last_search_nb_extern_pts;

  public:

    /// Constructor.
    /**
       \param p Parameters -- \b IN.
    */
    Cache_Search ( NOMAD::Parameters & p ) :
      NOMAD::Search              ( p , NOMAD::CACHE_SEARCH ) ,
      _last_search_nb_extern_pts ( 0                       )  {}
    
    /// Destructor.
    virtual ~Cache_Search ( void ) {}

    /// Reset.
    virtual void reset ( void ) { _last_search_nb_extern_pts = 0; }

    /// The cache search.
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
