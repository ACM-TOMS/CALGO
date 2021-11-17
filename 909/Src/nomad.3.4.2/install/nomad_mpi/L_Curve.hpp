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
  \file   L_Curve.hpp
  \brief  L_CURVE_TARGET stopping criterion (headers)
  \author Sebastien Le Digabel
  \date   2010-04-09
  \see    L_Curve.cpp
*/
#ifndef __L_CURVE__
#define __L_CURVE__

#include "Double.hpp"

namespace NOMAD {

  /// Class implementing the L_CURVE_TARGET stopping criterion.
  class L_Curve : private NOMAD::Uncopyable {

  private:
	  
    NOMAD::Double              _target;  ///< L_CURVE_TARGET parameter value.
    std::vector<NOMAD::Double> _f;       ///< List of objective values.
    std::vector<int          > _bbe;     ///< List of numbers of evaluations.

  public:

    /// Constructor.
    /**
       \param target L_CURVE_TARGET parameter value -- \b IN.
    */
    L_Curve ( const NOMAD::Double & target ) : _target ( target ) {}

    /// Destructor.
    virtual ~L_Curve ( void ) {}

    /// Insertion of a pair \c bbe/f in the lists \c _f and \c _bbe.
    /**
       \param bbe A new number of evaluations -- \b IN.
       \param f   A new objective value       -- \b IN.
    */
    void insert ( int bbe , const NOMAD::Double & f );

    /// Check the L_CURVE_TARGET stopping criterion.
    /**
       \param bbe An integer indicating a number of blackbox evaluations
                  -- \b IN.
       \return A boolean equal to \c true if the method detects that
               the target will not be reached after bbe evaluations.
    */
    bool check_stop ( int bbe ) const;
  };
}
#endif
