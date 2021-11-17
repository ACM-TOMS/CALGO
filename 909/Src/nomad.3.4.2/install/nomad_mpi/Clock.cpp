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
  \file   Clock.cpp
  \brief  Clock class (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-02
  \see    Clock.hpp
*/
#include "Clock.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
const double NOMAD::Clock::_D_CLOCKS_PER_SEC = static_cast<double>(CLOCKS_PER_SEC);

/*---------------------------------------------------------*/
/*  compute the wall-clock time (real time) elapsed since  */
/*  the construction of the Clock object                   */
/*---------------------------------------------------------*/
int NOMAD::Clock::get_real_time ( void ) const
{
  time_t t2;
  time  (&t2);
  return static_cast<int> (difftime ( t2 , _real_t0 ) );
}
