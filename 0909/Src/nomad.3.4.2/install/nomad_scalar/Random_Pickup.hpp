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
  \file   Random_Pickup.hpp
  \brief  Class for randomly pick up integers (headers)
  \author Sebastien Le Digabel
  \date   2010-04-07
  \see    Random_Pickup.cpp
*/
#ifndef __RANDOM_PICKUP__
#define __RANDOM_PICKUP__

#include <iostream>
#include <cstdlib>
#include "Uncopyable.hpp"

namespace NOMAD {
  
  /// Class for randomly pick up integers.
  /**
     - The integers are chosen in [0;n-1] and are distinct.
     - Example displaying 5 different integers in [0;4]:
     \code
     NOMAD::Random_Pickup rp(5);
     for ( int i = 0 ; i < 5 ; ++i )
       std::cout << rp.pickup() << std::endl;
     \endcode
  */
  class Random_Pickup : private NOMAD::Uncopyable {

  private:

    int   _n0;    ///< Initial value of \c n.
    int   _n;     ///< Current value of \c n.
    int * _elts;  ///< Elements that have not been chosen yet.

  public:

    /// Constructor.
    /**
       \param n -- The integer \c n defining the range
                   of values that can be picked up -- \b IN.
    */
    Random_Pickup ( int n );

    /// Destructor.
    virtual ~Random_Pickup ( void ) { delete [] _elts; }

    /// Reset.
    void reset ( void );

    /// Randomly pick up an element in [0;n-1].
    /**
       \return The element.
    */
    int pickup ( void );
    
    /// Cancel the last pick up.
    void cancel_last_pickup ( void );
    
  };
}

#endif
