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
  \file   Direction.hpp
  \brief  Polling direction (headers)
  \author Sebastien Le Digabel
  \date   2010-04-05
  \see    Direction.cpp
*/
#ifndef __DIRECTION__
#define __DIRECTION__

#include "Point.hpp"

namespace NOMAD {

  /// Class describing a polling direction.
  class Direction : public NOMAD::Point {

  private:

#ifdef DEBUG
    static int _cardinality;       ///< Number of NOMAD::Direction objects in memory.
    static int _max_cardinality;   ///< Max number of NOMAD::Direction objects in memory.
#endif

    NOMAD::direction_type _type;   ///< Type of direction.
    mutable int           _index;  ///< Direction index (used only for display).

  public:

#ifdef DEBUG
    /// Access to the number of NOMAD::Direction objects in memory.
    /**
       \return The number of NOMAD::Direction objects in memory.
    */
    static int get_cardinality ( void ) { return Direction::_cardinality; }

    /// Access to the max number of NOMAD::Direction objects in memory.
    /**
       \return The max number of NOMAD::Direction objects in memory.
    */
    static int get_max_cardinality ( void ) { return Direction::_max_cardinality; }
#endif

    /// Constructor #1.
    Direction ( void );

    /// Constructor #2.
    /**
       \param n Dimension -- \b IN.
       \param v Initial value for all coordinates -- \b IN.
       \param type Type of direction -- \b IN.
    */
    Direction ( int n , const NOMAD::Double & v , NOMAD::direction_type type );

    /// Constructor #3.
    /**
       \param x    Coordinates -- \b IN.
       \param type Type of direction -- \b IN.
    */
    Direction ( const NOMAD::Point & x , NOMAD::direction_type type );

    /// Copy constructor.
    /**
       \param d The copied object.
    */
    Direction ( const Direction & d );

    /// Destructor.
    virtual ~Direction ( void );

    /// Affectation operator.
    /**
       \param d The right-hand side object -- \b IN.
       \return \c *this as the result of the affectation.
    */
    Direction & operator = ( const Direction & d );

    /// Check if the direction is defined.
    /**
       \return A boolean equal to \c true if the direction has no defined type.
    */
    virtual bool is_defined ( void ) const { return _type != NOMAD::UNDEFINED_DIRECTION; }

    /// Clear the direction.
    virtual void clear ( void );

   /// Size of the direction in memory.
    /**
       \return An integer for the size of the direction in bytes.
    */    virtual int size_of ( void ) const
    {
      return NOMAD::Point::size_of() + sizeof(_type) + sizeof(_index);
    }

    /// Negation.
    /**
       The current object \c *this is not modified.
       \return A new direction equal to \c -*this.
    */
    const Direction operator - ( void ) const;
 
    /*---------------*/
    /*  GET methods  */
    /*---------------*/

    /// Access to the type of direction.
    /**
       \return The type of the direction.
    */
    NOMAD::direction_type get_type  ( void ) const { return _type; }
    
    /// Access to the direction index.
    /**
       \return The direction index.
    */
    int get_index ( void ) const { return _index; }

    /// Check if the direction is a MADS direction.
    /**
       \return A boolean equal to \c true if the direction is a MADS direction.
    */
    bool is_mads ( void ) const { return NOMAD::dir_is_mads ( _type ); }

   /// Check if the direction is a GPS direction.
    /**
       \return A boolean equal to \c true if the direction is a GPS direction.
    */
    bool is_gps ( void ) const { return NOMAD::dir_is_gps  ( _type ); }

    /*---------------*/
    /*  SET methods  */
    /*---------------*/

    /// Set the direction index.
    /**
       \param i The direction index -- \b IN.
    */
    void set_index ( int i ) const { _index = i; }

    /// Set the direction type.
    /**
       \param t The direction type -- \b IN.
    */
    void set_type ( NOMAD::direction_type t ) { _type  = t; }
    
    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
       \param sep A string that is used as a separator between the coordinates
                  -- \b IN.
       \param w   An integer indicating a width for the display of
                  each coordinate -- \b IN.
       \param lim Max number of coordinates to display -- \b IN.
    */
    virtual void display ( const NOMAD::Display & out ,
			   const std::string    & sep ,
			   int                    w   ,
			   int                    lim   ) const;
  };
  
  /// Display a NOMAD::Direction object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param d   The NOMAD::Direction object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display   & out ,
					      const NOMAD::Direction & d     )
  {
    d.display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
    return out;
  }
}

#endif
