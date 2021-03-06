/*

   MinDistance.
   Software package with several fast scalar, vector, and parallel
   implementations for computing the minimum distance of a random linear code.

   Copyright (C) 2017  Fernando Hernando (carrillf@mat.uji.es)
   Copyright (C) 2017  Francisco Igual (figual@ucm.es)
   Copyright (C) 2017  Gregorio Quintana (gquintan@uji.es)

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

*/


#include <stdio.h>
#include <stdlib.h>
#include "alg_common.h"


// ============================================================================
int get_next_combination( int * vecCombi,
    int numTotalElems, int numElemsToTake ) {
// Get new combination.
  int  diff, i, j, processedAll;

  if( vecCombi[ numElemsToTake - 1 ] < numTotalElems - 1 ) {
    // Most frequent case.
    vecCombi[ numElemsToTake - 1 ]++;
    processedAll = 0;
  } else {
    // General case.
    diff = numTotalElems - numElemsToTake;
    i = numElemsToTake;
    do {
      i--;
      vecCombi[ i ]++;
    } while( ( i > 0 )&&( vecCombi[ i ] > ( i + diff ) ) );
    for( j = i + 1; j < numElemsToTake; j++ ) {
      vecCombi[ j ] = vecCombi[ j - 1 ] + 1;
    }
    processedAll = ( vecCombi[ 0 ] > numTotalElems - numElemsToTake );
  }
  //// print_int_vector( "i_combi", vecCombi, numElemsToTake );
  return( processedAll );
}

