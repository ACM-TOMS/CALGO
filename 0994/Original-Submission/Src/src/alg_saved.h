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


#include "data_defs.h"

int generate_with_saved_alg(
    int cfg_num_saved_generators, int cfg_num_cores,
    GammaMatrix g, SavedCombi * vecSavedCombi,
    int numTotalElems, int numElemsToTake,
    int printMethod );

int process_prefix_with_saved_alg(
    int cfg_num_saved_generators, int cfg_num_cores,
    GammaMatrix g, SavedCombi * vecSavedCombi,
    int numTotalElems, int numElemsToTake,
    int lastValueInPrefix, uint32_t * vecAdditionPrefix,
    int printMethod );

void fill_structures_for_saved_data_for_1_s(
    int cfg_num_saved_generators,
    GammaMatrix g, SavedCombi * vecSavedCombi,
    int numTotalElems, int numElemsToTake,
    int printMethod );

