/****************************************************************************
 * RealPaver v. 0.4                                                         *
 *--------------------------------------------------------------------------*
 * Author: Laurent Granvilliers                                             *
 * Copyright (c) 1999-2003 Institut de Recherche en Informatique de Nantes  *
 * Copyright (c) 2004      Laboratoire d'Informatique de Nantes Atlantique  *
 *--------------------------------------------------------------------------*
 * RealPaver is distributed WITHOUT ANY WARRANTY. Read the associated       *
 * COPYRIGHT file for more details.                                         *
 *--------------------------------------------------------------------------*
 * interval.c                                                               *
 ****************************************************************************/


#include "interval.h"

/* interval constants to be initialized */
IBItv IBConstantPi;
IBItv IBConstant_2_Pi;
IBItv IBConstant_Pi_2;


void IBInitIntervalConstants() {
/***************************************************************************
*  
*/
  IBSetToPiI(IBConstantPi);
  IBMulRposIinternal(IBConstant_2_Pi,2.0,IBConstantPi);
  IBMulRposIinternal(IBConstant_Pi_2,0.5,IBConstantPi);
}
