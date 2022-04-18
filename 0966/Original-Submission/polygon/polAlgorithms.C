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


#include "polAlgorithms.h"

/**
 * Computes a "pseudo" angle of an edge from 'orig' to 'dest'. Consider two
 * angles 'a' and 'b': 
 *  - if 'a' > 'b' <----> pseudo('a') > pseudo('b')  
 *  - if 'a' < 'b' <----> pseudo('a') < pseudo('b')  
 *  - if 'a' = 'b' <----> pseudo('a') = pseudo('b')  
 */
RT getPseudoAngle (Point orig, Point dest) {
    Line line(Point(0,0), (Segment(orig, dest)).direction());

    if (dest.y() > orig.y()) {
        if (dest.x() >= orig.x()) {
            // the angle is at the first quadrant
            RT xCoord = line.x_at_y(RT(1));
            if (xCoord <= RT(1)) {
                return RT(2) - xCoord;
            } else {
                return RT(0) + line.y_at_x(RT(1));
            }
        } else {
            // the angle is at the second quadrant
            RT xCoord = line.x_at_y(RT(1));
            if (xCoord >= RT(-1)) {
                return RT(2) - xCoord;
            } else {
                return RT(4) - line.y_at_x(RT(-1));
            }
        }
    } else if (dest.y() < orig.y()) {
        if (dest.x() >= orig.x()) {
            // the angle is at the fourth quadrant
            RT xCoord = line.x_at_y(RT(-1));
            if (xCoord <= RT(1)) {
                return RT(6) + xCoord;
            } else {
                return RT(8) + line.y_at_x(RT(1));
            }
        } else {
            // the angle is at the third quadrant
            RT xCoord = line.x_at_y(RT(-1));
            if (xCoord >= RT(-1)) {
                return RT(6) + xCoord;
            } else {
                return RT(4) - line.y_at_x(RT(-1));
            }
        }
    } else {
        if (dest.x() >= orig.x()) {
            return RT(0);
        } else {
            return RT(4);
        }
    }
}

/**
 * Computes de angle between two lines defined by points 'a', 'b' and 'c'.
 */
RT getAngle(Point a, Point b, Point c) {
    Segment s1(b, a);
    Segment s2(b, c);

    RT angle = getPseudoAngle(b, a) - getPseudoAngle(b, c);

    if (angle < RT(0)) {
        angle = RT(0) - angle;
    }

    if (angle > RT(4)) {
        angle = RT(8) - angle;
    }

    return angle;
}

/**
 * Finds Alpha angle.
 */
RT getAlpha(RT alphaAnt, Point a, Point b, Point c) {
    if (CGAL::collinear(b, a, c)) {
        return alphaAnt;
    }
    RT angle = getAngle(a, b, c);

    if (CGAL::left_turn(b, a, c)) {
        return alphaAnt + angle;
    } else {
        // it is CGAL::right_turn(a, b, c)
        return alphaAnt - angle;
    }
}

/**
 * Simulates an point in the infinite direction start -> dir getting the
 * intersection of this ray with the bounding box.
 */
Point pointOnInfDir(Point start, Point dir, Polygon pbox) {
    Ray b = Ray(start, dir);

    Polygon::Edge_const_iterator eIt = pbox.edges_begin();
    for (; eIt != pbox.edges_end(); ++eIt) {
        CGAL::Object obj = CGAL::intersection(*eIt, b);
        if (const Point *point = CGAL::object_cast<Point>(&obj)) {
            return *point;
        }
    }
    // shouldn't get here. but if so, return an valid point anyway
    return start;
}

/**
 * Cleans a polygon, removing equal and colinear points.
 */
PolygonExt cleanPol(PolygonExt pol) {
    // remove equal points
    PolygonExt newPol;
    Polygon::Vertex_iterator vIt = pol.vertices_begin();
    Point oldVertex = (*vIt);
    Point firstVertex = (*vIt);
    newPol.push_back(oldVertex);
    vIt++;
    for (; vIt != pol.vertices_end(); vIt++) {
        if (oldVertex != (*vIt) && firstVertex != (*vIt)) {
            oldVertex = (*vIt);
            newPol.push_back(oldVertex);
        }
    }
    // remove colinear vertices
    PolygonExt otherPol;
    Polygon::Vertex_circulator vIminus1 = newPol.vertices_circulator();
    Polygon::Vertex_iterator vI = newPol.vertices_begin();
    Polygon::Vertex_circulator vIplus1 = newPol.vertices_circulator();
    for(;(*vIplus1) != (*vI); vIplus1++);
    for(;(*vIminus1) != (*vI); vIminus1++);
    vIplus1++;
    vIminus1--;
    for (; vI != newPol.vertices_end(); vI++) {
        if (!CGAL::collinear(*vIminus1, *vI, *vIplus1)) {
            otherPol.push_back(*vI);
        }
        vIplus1++;
        vIminus1++;
    }
    return otherPol;
}
