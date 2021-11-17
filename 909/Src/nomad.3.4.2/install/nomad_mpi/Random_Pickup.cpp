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
  \file   Random_Pickup.cpp
  \brief  Class for randomly pick up integers (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-07
  \see    Random_Pickup.hpp
*/
#include "Random_Pickup.hpp"

/*---------------------------------------------------------*/
/*                         constructor                     */
/*---------------------------------------------------------*/
NOMAD::Random_Pickup::Random_Pickup ( int n )
  : _n0   ( n          ) ,
    _n    ( n          ) ,
    _elts ( new int[n] )
{
  for ( int i = 0 ; i < n ; ++i )
    _elts[i] = i;
}

/*---------------------------------------------------------*/
/*                           reset                         */
/*---------------------------------------------------------*/
void NOMAD::Random_Pickup::reset ( void )
{
  _n = _n0;
  for ( int i = 0 ; i < _n ; ++i )
    _elts[i] = i;
}

/*---------------------------------------------------------*/
/*                randomly pick up an element              */
/*---------------------------------------------------------*/
int NOMAD::Random_Pickup::pickup ( void )
{
  if ( _n == 0 )
    return 0;
  int ind = rand()%_n;
  int tmp = _elts[ind];
  if ( ind < _n - 1 ) {
    _elts[ind ] = _elts[_n-1];
    _elts[_n-1] = tmp;
  }
  --_n;
  return tmp;
}

/*---------------------------------------------------------*/
/*                   cancel the last pick up               */
/*---------------------------------------------------------*/
void NOMAD::Random_Pickup::cancel_last_pickup ( void )
{
  if ( _n < _n0 )
    ++_n;
}
