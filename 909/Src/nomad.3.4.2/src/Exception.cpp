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
  \file   Exception.cpp
  \brief  custom class for exceptions (implementation)
  \author Sebastien Le Digabel
  \date   2010-03-29
  \see    Exception.hpp
*/
#include "Exception.hpp"

/*----------------------------------------------------------------*/
/*                     NOMAD::Exception::what()                   */
/*----------------------------------------------------------------*/
const char * NOMAD::Exception::what ( void ) const throw()
{
  std::ostringstream oss;
  oss << "NOMAD::Exception thrown (" << _file << ", " << _line << ")";
  if ( !_what.empty() )
    oss << " " << _what;
  _what = oss.str();
  return _what.c_str();
}
