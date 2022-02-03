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
  \file   VNS_Search.hpp
  \brief  VNS search (headers)
  \author Sebastien Le Digabel
  \date   2010-04-12
  \see    VNS_Search.cpp
*/
#ifndef __VNS_SEARCH__
#define __VNS_SEARCH__

#include "Mads.hpp"
#include "Search.hpp"

namespace NOMAD {

  /// Variable Neighborhood Search (VNS) search.
  class VNS_Search : public NOMAD::Search , private NOMAD::Uncopyable {

  private:

    int                           _k;  ///< VNS neighborhood parameter.
    int                       _k_max;  ///< Maximum value of \c _k.
    int                _halton_index;  ///< Halton index used for shaking directions.
    const NOMAD::Eval_Point * _old_x;  ///< Previous reference point (updates \c _k).
    
    /// Search frequency.
    /**
       \c _nb_performed[ell] corresponds to the number of times that
       the VNS search has been performed on a mesh of a size \c ell.
    */
    int _nb_performed [NOMAD::L_LIMITS+1];

  public:

    /// Constructor.
    /**
       \param p Parameters -- \b IN.
    */
    VNS_Search ( NOMAD::Parameters & p )
      : NOMAD::Search ( p , NOMAD::VNS_SEARCH     ) ,
	_k            ( 1                         ) ,
	_k_max        ( 1                         ) ,
	_halton_index ( NOMAD::VNS_HALTON_INDEX_0 ) ,
	_old_x        ( NULL                      )   {}
    
    /// Destructor.
    virtual ~VNS_Search ( void ) {}

    /// Reset.
    virtual void reset ( void );

    /// The VNS search.
    /**
       Principle: x --[shaking(k)]--> x' --[descent]--> x" .
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
  
    /// Access to the search frequency.
    /**
       \param ell The mesh index -- \b IN.
       \return Number of times that the VNS search has been performed on a given mesh.
    */
    int get_nb_performed ( int ell ) const
    {
      return ( ell < 0 ) ? 0 : _nb_performed[ell];
    }
  };
}

#endif
