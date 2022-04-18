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
  \file   Cache_Point.hpp
  \brief  Point stored in the cache (headers)
  \author Sebastien Le Digabel
  \date   2010-04-08
  \see    Cache_Point.cpp
*/
#ifndef __CACHE_POINT__
#define __CACHE_POINT__

#include "Eval_Point.hpp"

namespace NOMAD {

  /// Class for the representation of NOMAD::Eval_Point objects stored in the cache.
  class Cache_Point : public Set_Element<NOMAD::Eval_Point> {

  private:

    /// Affectation operator.
    /**
       \param cp The right-hand side object -- \b IN.
       \return \c *this as the result of the affectation.
    */
    Cache_Point & operator = ( const Cache_Point & cp );

  public:

    /// Constructor.
    /**
       \param x A pointer to the NOMAD::Eval_Point object
                that is stored in the cache -- \b IN.
    */
    explicit Cache_Point ( const NOMAD::Eval_Point * x ) :
      Set_Element<NOMAD::Eval_Point> ( x ) {}

    /// Copy constructor.
    /**
       \param cp The copied object -- \b IN.
    */
    Cache_Point ( const Cache_Point & cp )
      : Set_Element<NOMAD::Eval_Point> ( cp ) {}

    /// Destructor.
    virtual ~Cache_Point ( void ) {}

    /// Comparison operator.
    /**
       \param  cp The right-hand side object.
       \return A boolean equal to \c true if \c *this \c < \c cp.
    */
    virtual bool operator < ( const Set_Element<NOMAD::Eval_Point> & cp ) const;

    /// Access to the point.
    /**
       \return A pointer to the NOMAD::Eval_Point stored in the cache.
    */
    const NOMAD::Eval_Point * get_point ( void ) const { return get_element(); }
  };
}

#endif
