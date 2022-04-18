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
  \file   Clock.hpp
  \brief  Clock class (headers)
  \author Sebastien Le Digabel
  \date   2010-04-02
  \see    Clock.cpp
*/
#ifndef __CLOCK__
#define __CLOCK__

#include <ctime>

namespace NOMAD {

  /// Clock class.
  /**
     Time measurement.\n\n
     \b Example:
     \code
     Clock c;

     // some instructions here

     std::cout << "elapsed real time = " << c.get_real_time() << std::endl;
     std::cout << "elapsed cpu time  = " << c.get_cpu_time()  << std::endl;
     \endcode
  */
  class Clock {

  private:

    time_t              _real_t0;          ///< Wall clock time measurement.
    clock_t             _cpu_t0;           ///< CPU time measurement.
    static const double _D_CLOCKS_PER_SEC; ///< System constant for CPU time measurement.

  public:

    /// Constructor.
    Clock ( void ) : _cpu_t0 ( clock() ) { time (&_real_t0); }

    /// Copy constructor.
    /**
       \param c The copied object -- \b IN.
    */
    Clock ( const Clock & c ) : _real_t0 ( c._real_t0 ) , _cpu_t0 ( c._cpu_t0 ) {}

    /// Affectation operator.
    /**
       \param  c The right-hand side object -- \b IN.
       \return \c *this as the result of the affectation.
    */
    Clock & operator = ( const Clock & c )
    {
      _real_t0 = c._real_t0;
      _cpu_t0  = c._cpu_t0;
      return *this;
    }

    /// Destructor.
    virtual ~Clock ( void ) {}

    /// Reset the clock.
    void reset ( void )
    {
      time ( &_real_t0 );
      _cpu_t0 = clock();
    }

    /// Get wall clock time.
    /**
       \return The wall clock time.
    */
    int get_real_time ( void ) const;

    /// Get the cpu time.
    /**
       \return The cpu time.
    */
    double get_cpu_time ( void ) const
    {
      return ( clock() - _cpu_t0 ) / _D_CLOCKS_PER_SEC;
    }
  };
}

#endif
