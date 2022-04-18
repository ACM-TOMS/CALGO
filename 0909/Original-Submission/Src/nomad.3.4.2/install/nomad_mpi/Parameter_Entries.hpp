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
  \file   Parameter_Entries.hpp
  \brief  Parameter entries (headers)
  \author Sebastien Le Digabel
  \date   2010-04-05
  \see    Parameter_Entries.cpp
*/
#ifndef __PARAMETER_ENTRIES__
#define __PARAMETER_ENTRIES__

#include "Parameter_Entry.hpp"

namespace NOMAD {

  /// Parameter entries.
  /**
     - Objects of this class store NOMAD::Parameter_Entry objects.
     - One NOMAD::Parameter_Entries object summarizes an entire parameters file.
  */
  class Parameter_Entries : private NOMAD::Uncopyable {

  private:

    /// List of NOMAD::Parameter_Entry objects (the entries).
    std::multiset<NOMAD::Parameter_Entry *, NOMAD::Parameter_Entry_Comp> _entries;

  public:

    /// Constructor.
    explicit Parameter_Entries ( void ) {}

    /// Destructor.
    virtual ~Parameter_Entries ( void );

    /// Find a specific entry in a set.
    /**
       \param  name The name of the wanted NOMAD::Parameter_Entry object -- \b IN.
       \return      A pointer to the NOMAD::Parameter_Entry object if it
                    has been found in the list of entries,
                    or \c NULL otherwise.
    */
    NOMAD::Parameter_Entry * find ( const std::string & name ) const;
  
    /// Insert a new entry in the list of entries.
    /**
       \param entry A pointer to the new NOMAD::Parameter_Entry object -- \b IN.
    */
    void insert ( NOMAD::Parameter_Entry * entry );

    /// Find a non-interpreted entry.
    /**
       \return A pointer to the first NOMAD::Parameter_Entry that has not been
               interpreted so far,
               or \c NULL if all entries have already been interpreted.
    */
    NOMAD::Parameter_Entry * find_non_interpreted ( void ) const;

    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */    void display ( const NOMAD::Display & out ) const;
  };

  /// Display a NOMAD::Parameter_Entries object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param e   The NOMAD::Parameter_Entries object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display           & out ,
					      const NOMAD::Parameter_Entries & e     )
  {
    e.display ( out );
    return out;
  }
}

#endif
