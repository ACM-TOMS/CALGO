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


#ifndef POLYGON_EXT_H
#define POLYGON_EXT_H

#include <CGAL/Gmpq.h>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/intersections.h>
#include <list>
#include <math.h>

typedef CGAL::Gmpq RT;
typedef CGAL::Cartesian<RT> Extended_kernel;


typedef CGAL::Point_2<Extended_kernel> Point;
typedef CGAL::Vector_2<Extended_kernel> Vector;
typedef CGAL::Segment_2<Extended_kernel> Segment;
typedef CGAL::Ray_2<Extended_kernel> Ray;
typedef CGAL::Polygon_2<Extended_kernel, std::list<Point> > Polygon;
typedef CGAL::Polygon_with_holes_2<Extended_kernel, std::list<Point> > PolygonWithHoles;
typedef CGAL::Line_2<Extended_kernel> Line;

class PolygonExt : public Polygon{
    public:
       PolygonExt() : Polygon() {}
       PolygonExt(Polygon pol) : Polygon(pol) {}
       
       /** 
        * Finds the closest point on the boundary to z.
        */
       Point getClosestOnBoundary(Point);
       
       /**
        * Computes the visibility polygon of a point z in P.
        * Algorithm: B. Joe and R. B. Simpson. Visibility of a simple polygon from a point. 
        * Report CS-85-38, Dept. Math. Comput. Sci., Drexel Univ., Philadelphia, PA, 1985.
        */
       PolygonExt getVisibility(Point);
};

#endif

