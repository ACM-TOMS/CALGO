/*****************************************************************************
 * This code is part of Art Gallery Solver (AGSol) Package, which aims the 
 * resolution of the Art Gallery Problem With Point Guards.
 *
 * This software version (1.0.2) has been tested under and is compatible with 
 * CGAL 3.9 and GLPK 4.52.
 *
 * Authors:
 *  Davi Colli Tozoni - davi.tozoni@gmail.com
 *  Marcelo Castilho Couto - coutomarcelo@gmail.com
 *
 * AGSol Concept and Design: 
 *  Davi Colli Tozoni, Marcelo Castilho Couto, Pedro Jussieu de Rezende & Cid 
 * Carvalho de Souza.
 * 
 * Other information can be found at:
 *  http://www.ic.unicamp.br/~cid/Problem-instances/Art-Gallery/index.html
 *
 * --
 * 
 *  This program is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation, either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.
 *      
 ****************************************************************************/


#ifndef AVP_LIGHT_GRID_H
#define AVP_LIGHT_GRID_H

#include "IGrid.h"
#include "Arrangement.h"

class AVPLightGrid : public IGrid {
    private:
        struct rusage ruseTime;

        Arrangement* _arr;
        MyObserver* _obs;

        int _vertexNum;
        int _edgeNum;
        int _faceNum;

    public:
        /**
         * Destrutor. Frees arrangement's and observer's memory
         */
        ~AVPLightGrid();

        /**
         * Constructs the grid.
         */
        virtual void makeGrid();
        
        /**
         * Adds edges from new visibility polygons in the arrangement.
         */
        void addToGrid(std::vector<PolygonExt> vet);
};

#endif
