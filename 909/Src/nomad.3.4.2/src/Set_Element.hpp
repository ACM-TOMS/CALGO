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
  \file   Set_Element.hpp
  \brief  Element of a set (headers)
  \author Sebastien Le Digabel
  \date   2010-04-12
*/
#ifndef __SET_ELEMENT__
#define __SET_ELEMENT__

namespace NOMAD {

  // forward declarations:
  class Double;
  class Eval_Point;

  /// Generic class for elements of a \c std::set.
  /**
     This is an abstract class (it is not possible to create NOMAD::Set_Element objects).
  */
  template <class T>
  class Set_Element {

  private:

#ifdef DEBUG
    static int _cardinality;     ///< Number of NOMAD::Set_Element objects in memory.
    static int _max_cardinality; ///< Max number of NOMAd::Set_Element objects in memory.
#endif

    const T * _el; ///< A pointer to the element.

    /// Affectation operator.
    /**
       \param se The right-hand side object -- \b IN.
    */
    Set_Element & operator = ( const Set_Element & se );

  public:
  
    /// Constructor.
    /**
       \param el A pointer on the element -- \b IN.
    */
    explicit Set_Element ( const T * el ) : _el ( el ) {
#ifdef DEBUG
      ++Set_Element::_cardinality;
      if ( Set_Element::_cardinality > Set_Element::_max_cardinality )
	++Set_Element::_max_cardinality;
#endif
    }
  
    /// Copy constructor.
    /**
       \param sp The copied object -- \b IN.
    */
    explicit Set_Element ( const Set_Element & sp ): _el ( sp._el ) {
#ifdef DEBUG
      ++Set_Element::_cardinality;
      if ( Set_Element::_cardinality > Set_Element::_max_cardinality )
	++Set_Element::_max_cardinality;
#endif
    }

    /// Destructor.
    virtual ~Set_Element ( void ) {
#ifdef DEBUG
      --Set_Element::_cardinality;
#endif
    }

    /// Specific NOMAD::Priority_Eval_Point elements of comparison.
    /**
       - Only NOMAD::Priority_Eval_Point::get_priority_criteria() does something.
       - \see Priority_Eval_Point.hpp .
       \param c1 A real -- \b IN.
       \param c2 A real -- \b IN.
       \param c3 A real -- \b IN.
       \param c4 A real -- \b IN.
    */
    virtual void get_priority_criteria ( NOMAD::Double & c1 ,
					 NOMAD::Double & c2 ,
					 NOMAD::Double & c3 ,
					 NOMAD::Double & c4   ) const {}
    /// Comparison operator.
    /**
       - Has to be implemented by every NOMAD::Set_Element subclass.
       - Pure virtual method.
       \param se The right-hand side object -- \b IN.
    */
    virtual bool operator < ( const Set_Element & se ) const = 0;
    
    /// Access to the element.
    /**
       \return A pointer to the element.
    */
    const T * get_element ( void ) const { return _el; }
    
    /// Set an element.
    /**
       \param el A pointer to the element -- \b IN.
    */
    void set_element ( const T * el ) { _el = el; }
    
#ifdef DEBUG

    /// Access to the number of NOMAD::Set_Element objects in memory.
    /**
       \return Number of NOMAD::Set_Element objects in memory.
    */
    static int get_cardinality ( void ) { return Set_Element::_cardinality; }

    /// Access to the max number of NOMAD::Set_Element objects in memory.
    /**
       \return Max number of NOMAD::Set_Element objects in memory.
    */
    static int get_max_cardinality ( void ) { return Set_Element::_max_cardinality; }
    
#endif
  };

#ifdef DEBUG

  /// Initialization of _cardinality.
  template<class T> int NOMAD::Set_Element<T>::_cardinality = 0;

  /// Initialization of _max_cardinality.
  template<class T> int NOMAD::Set_Element<T>::_max_cardinality = 0;
#endif
}

#endif
