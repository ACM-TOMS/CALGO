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


#ifndef IGRID_H
#define IGRID_H

#include "PolygonExt.h"
#include "PolygonWithHolesExt.h"

#include <time.h>
#include <stdlib.h>
#include <malloc.h>

/* Procedure for measuring runtime. Implemented by Ionut Aaron (Brown,GSIA). */
#ifndef __SYSTEM_TIME_H
#define __SYSTEM_TIME_H

/* Ionut Aron, ia@cs.brown.edu (1999-2001)
 * created: August, 1999
 * last updated: December 20, 2001
 */

#include <sys/resource.h>
#include <sys/types.h>
#include <time.h>

extern int getrusage();

// this macro produces a time equal to the one produced by clock(),
// but does not suffer from the wraparound problem of clock()

#define CPUTIME(ruse) (\
  getrusage(RUSAGE_SELF,&ruse),\
  ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + \
  1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec)\
)

/* usage of the CPUTIME timing routine:
   ===================================

   struct rusage ruse;
   double t0 = CPUTIME(ruse);
   ...
   double t1 = CPUTIME(ruse);
   printf("CPU time: %f secs.\n", t1-t0);
*/

// global variable
// struct rusage ruseTime;

#endif

class IGrid {
    protected:
        PolygonWithHolesExt _polygon;
        std::vector<PolygonExt> _visPol;
        std::vector<Point> _grid;

        bool _visOk;
        bool _debug;
        int _gridMem;
        double _gridTime;

        std::vector<PolygonExt> _struct1;
        std::vector<PolygonExt> _struct2;

    public:
        /**
         * Constructor and Destructor.
         */
        IGrid() { _visOk = false; _debug = false; _gridMem = 0; }
        virtual ~IGrid() { }

        /**
         * Sets polygon to be treated.
         */
        void setPolygon(PolygonWithHolesExt polygon);
        
        /**
         * Copies the visibility polygons to a structure inside IGrid class.
         */
        void setVisibilityPolygons(std::vector<PolygonExt> visPol);

        /**
         * Constructs the grid.
         */
        virtual void makeGrid() = 0;
        
        /**
         * Returns memory used by the grid structure.
         */
        int getGridMem() { return _gridMem; }

        /**
         * Adds a point 'p' to the current grid.
         */
        void addGridPoint(Point p);

        /**
         * Returns grid points.
         */
        std::vector<Point> getGridPoints() { return _grid; }

        /**
         * Returns first struct.
         */
        std::vector<PolygonExt> getStruct1() { return _struct1; }

        /**
         * Returns second struct.
         */
        std::vector<PolygonExt> getStruct2() { return _struct2; }
        
};

#endif

