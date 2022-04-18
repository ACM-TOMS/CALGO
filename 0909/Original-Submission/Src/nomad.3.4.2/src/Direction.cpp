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
  \file   Direction.cpp
  \brief  Polling direction (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-05
  \see    Direction.hpp
*/
#include "Direction.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
#ifdef DEBUG
int NOMAD::Direction::_cardinality     = 0;
int NOMAD::Direction::_max_cardinality = 0;
#endif

/*---------------------------------------------------------*/
/*                       constructor 1                     */
/*---------------------------------------------------------*/
NOMAD::Direction::Direction ( void )
  : NOMAD::Point  (                            ) ,
    _type         ( NOMAD::UNDEFINED_DIRECTION ) ,
    _index        ( -1                         )
{
#ifdef DEBUG
  ++NOMAD::Direction::_cardinality;
  if ( NOMAD::Direction::_cardinality > NOMAD::Direction::_max_cardinality )
    ++NOMAD::Direction::_max_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                       constructor 2                     */
/*---------------------------------------------------------*/
NOMAD::Direction::Direction ( int                     n    ,
			      const NOMAD::Double   & v    ,
			      NOMAD::direction_type   type   )
  : NOMAD::Point  ( n , v ) ,
    _type         ( type  ) ,
    _index        ( -1    )
{
#ifdef DEBUG
  ++NOMAD::Direction::_cardinality;
  if ( NOMAD::Direction::_cardinality > NOMAD::Direction::_max_cardinality )
    ++NOMAD::Direction::_max_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                       constructor 3                     */
/*---------------------------------------------------------*/
NOMAD::Direction::Direction ( const NOMAD::Point    & x    , 
			      NOMAD::direction_type   type   )
  : NOMAD::Point  ( x    ) ,
    _type         ( type ) ,
    _index        ( -1   )
{
#ifdef DEBUG
  ++NOMAD::Direction::_cardinality;
  if ( NOMAD::Direction::_cardinality > NOMAD::Direction::_max_cardinality )
    ++NOMAD::Direction::_max_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                      copy constructor                   */
/*---------------------------------------------------------*/
NOMAD::Direction::Direction ( const Direction & d )
  : NOMAD::Point  ( d        ) ,
    _type         ( d._type  ) ,
    _index        ( d._index )
{
#ifdef DEBUG
  ++NOMAD::Direction::_cardinality;
  if ( NOMAD::Direction::_cardinality > NOMAD::Direction::_max_cardinality )
    ++NOMAD::Direction::_max_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                         destructor                      */
/*---------------------------------------------------------*/
NOMAD::Direction::~Direction ( void )
{
#ifdef DEBUG
  --NOMAD::Direction::_cardinality;
#endif
}

/*---------------------------------------------------------*/
/*                  affectation operator                   */
/*---------------------------------------------------------*/
NOMAD::Direction & NOMAD::Direction::operator = ( const NOMAD::Direction & d )
{
  if ( this == &d )
    return *this;

  NOMAD::Point::operator = ( d );
  
  _type  = d._type;
  _index = d._index;

  return *this;
}

/*---------------------------------------------------------*/
/*                            clear                        */
/*---------------------------------------------------------*/
void NOMAD::Direction::clear ( void )
{
  NOMAD::Point::clear();
  _type  = NOMAD::UNDEFINED_DIRECTION;
  _index = -1;
}

/*-----------------------------------------------------------*/
/*                           negation                        */
/*-----------------------------------------------------------*/
const NOMAD::Direction NOMAD::Direction::operator - ( void ) const
{
  return NOMAD::Direction ( this->NOMAD::Point::operator-() , _type );
}

/*---------------------------------------------------------*/
/*                          display                        */
/*---------------------------------------------------------*/
void NOMAD::Direction::display ( const NOMAD::Display & out ,
				 const std::string    & sep ,
				 int                    w   ,
				 int                    lim   ) const
{
  if ( is_defined() ) {
    out << "( ";
    NOMAD::Point::display ( out , sep , w , lim );
    out << " ) " << _type;
  }
  else
    out << "undefined";
}
