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


#ifndef POL_ALGO_H
#define POL_ALGO_H

#include "PolygonExt.h"

/**
 * Computes a "pseudo" angle of an edge from 'orig' to 'dest'. Consider two
 * angles 'a' and 'b': 
 *  - if 'a' > 'b' <----> pseudo('a') > pseudo('b')  
 *  - if 'a' < 'b' <----> pseudo('a') < pseudo('b')  
 *  - if 'a' = 'b' <----> pseudo('a') = pseudo('b')  
 */
RT getPseudoAngle (Point orig, Point dest);

/**
 * Computes de angle between two lines defined by points 'a', 'b' and 'c'.
 */
RT getAngle(Point a, Point b, Point c);

/**
 * Finds Alpha angle.
 */
RT getAlpha(RT alphaAnt, Point a, Point b, Point c);

/**
 * Simulates an point in the infinite direction start -> dir getting the
 * intersection of this ray with the bounding box.
 */
Point pointOnInfDir(Point start, Point dir, Polygon pbox);

/**
 * Cleans a polygon, removing equal and colinear points.
 */
PolygonExt cleanPol(PolygonExt pol);

#endif
