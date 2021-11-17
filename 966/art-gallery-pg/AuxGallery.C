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


#include "AuxGallery.h"

/**
 * If a hole is a non simple polygon, split the same in simple ones.
 */
vector<PolygonExt> splitNonSimple(PolygonExt pol) {
    vector<PolygonExt> vet;
    Vertex_iterator init2 = pol.vertices_begin();
    Vertex_iterator vit2;
    bool found = false;
    
    PolygonExt pol1, pol2, pol2tmp;
    for(Vertex_iterator vit = pol.vertices_begin(); vit != pol.vertices_end(); ++vit){
        init2 = vit;
        init2++;

        for(vit2 = init2; vit2 != pol.vertices_end(); ++vit2){
            pol2.push_back(*vit2);
            if ((*vit2) == (*vit)) {
                found = true;
                break;
            }
        }

        if (found == true) {
            break;
        }
        else {
            pol2.clear();
        }
        
        pol1.push_back(*vit);
    }
    
    for(Vertex_iterator vit = vit2; vit != pol.vertices_end(); ++vit){
        pol1.push_back(*vit);
    }

    vet.push_back(pol1);
    vet.push_back(pol2);
    
    return vet;
}

/**
 * Gets file name, whithout the whole path.
 */
string getFileName(string polFileName){
    std::string str = polFileName;
    int found;

    found = str.find_last_of(".");
    
    return str.substr(0, found);
}

/**
 * Finds an internal point of a Polygon With Holes.
 */
Point internalPoint(PolygonWithHoles polWHoles) {
    Polygon pol = polWHoles.outer_boundary();
    RT xCentral = 0;
    RT yCentral = 0;
    RT numVertex = RT(pol.size());

    /* First thing is to identify a convex vertex v and its adjacents. */
    Polygon::Vertex_circulator vIminus1 = pol.vertices_circulator();
    Polygon::Vertex_circulator vI = pol.vertices_circulator();
    Polygon::Vertex_circulator vIplus1 = pol.vertices_circulator();
    for(;(*vIplus1) != (*vI); vIplus1++);
    for(;(*vIminus1) != (*vI); vIminus1--);
    vIplus1++;
    vIminus1--;
    Point convexVertex = *vI;

    do {
        if (CGAL::left_turn(*vIminus1, *vI, *vIplus1)) {
            convexVertex = *vI;
            break;
        }
        vIplus1++;
        vIminus1++;
        vI++;
    } while (*vI != convexVertex);

    Line lineAB(*vIplus1, *vIminus1);
    Polygon polAVB;
    polAVB.push_back(*vIminus1);
    polAVB.push_back(*vI);
    polAVB.push_back(*vIplus1);

    RT maxDist = RT(-1);
    Point maxVertex;

    /* For each other vertex q, if q is inside polAVB, computes distance to q
     * (orthogonal to ab).  If distance is a new maximum, save the point. */
    vI++;
    vI++;
    for (; *vI != *vIminus1; vI++) {
        if (polAVB.bounded_side(*vI) == CGAL::ON_BOUNDED_SIDE) {
            RT dist = CGAL::squared_distance(*vI, lineAB);
            if (dist > maxDist) {
                maxVertex = *vI;
                maxDist = dist;
            }
        }
    }

    /* Verify holes */    
    for(PolygonWithHoles::Hole_const_iterator holeIt = polWHoles.holes_begin(); holeIt != polWHoles.holes_end(); ++holeIt) {
        for (Polygon::Vertex_iterator vIt = (*holeIt).vertices_begin(); vIt != (*holeIt).vertices_end(); vIt++) {
            if (polAVB.bounded_side(*vIt) == CGAL::ON_BOUNDED_SIDE) {
                RT dist = CGAL::squared_distance(*vIt, lineAB);
                if (dist > maxDist) {
                    maxVertex = *vIt;
                    maxDist = dist;
                }
            }
        }
    }

    /* If no point is inside, return midpoint of ab, or centroid of polAVP.
       Otherwise, if some point inside qv is internal, return its midpoint. */
    if (maxDist == RT(-1)) {
        xCentral = ((convexVertex).x() + (*vIminus1).x() + (*vIplus1).x()) / RT(3);
        yCentral = ((convexVertex).y() + (*vIminus1).y() + (*vIplus1).y()) / RT(3);
    } else {
        xCentral = ((convexVertex).x() + (maxVertex).x()) / RT(2);
        yCentral = ((convexVertex).y() + (maxVertex).y()) / RT(2);
    }

    return Point(xCentral,yCentral);
}

/**
 * Finds an internal point of a Polygon Without Holes.
 */
Point internalPoint(Polygon pol){
    RT xCentral = 0;
    RT yCentral = 0;
    RT numVertex = RT(pol.size());

    /* First thing is to identify a convex vertex v and its adjacents. */
    Polygon::Vertex_circulator vIminus1 = pol.vertices_circulator();
    Polygon::Vertex_circulator vI = pol.vertices_circulator();
    Polygon::Vertex_circulator vIplus1 = pol.vertices_circulator();
    for(;(*vIplus1) != (*vI); vIplus1++);
    for(;(*vIminus1) != (*vI); vIminus1--);
    vIplus1++;
    vIminus1--;
    Point convexVertex = *vI;

    do {
        if (CGAL::left_turn(*vIminus1, *vI, *vIplus1)) {
            convexVertex = *vI;
            break;
        }
        vIplus1++;
        vIminus1++;
        vI++;
    } while (*vI != convexVertex);

    Line lineAB(*vIplus1, *vIminus1);
    Polygon polAVB;
    polAVB.push_back(*vIminus1);
    polAVB.push_back(*vI);
    polAVB.push_back(*vIplus1);

    RT maxDist = RT(-1);
    Point maxVertex;

    /* For each other vertex q, if q is inside polAVB, computes distance to q
     * (orthogonal to ab).  If distance is a new maximum, save the point. */
    vI++;
    vI++;
    for (; *vI != *vIminus1; vI++) {
        if (polAVB.bounded_side(*vI) == CGAL::ON_BOUNDED_SIDE) {
            RT dist = CGAL::squared_distance(*vI, lineAB);
                if (dist > maxDist) {
                    maxVertex = *vI;
                    maxDist = dist;
                }
        }
    }

    /* If no point is inside, return midpoint of ab, or centroid of polAVP.
       Otherwise, if some point inside qv is internal, return its midpoint. */
    if (maxDist == RT(-1)) {
        xCentral = ((convexVertex).x() + (*vIminus1).x() + (*vIplus1).x()) / RT(3);
        yCentral = ((convexVertex).y() + (*vIminus1).y() + (*vIplus1).y()) / RT(3);
    } else {
        xCentral = ((convexVertex).x() + (maxVertex).x()) / RT(2);
        yCentral = ((convexVertex).y() + (maxVertex).y()) / RT(2);
    }
    return Point(xCentral,yCentral);

}

/**
 * Finds the size of the fraction x.
 */
int fractionSize(RT x) {
    std::string s;
    std::stringstream out;
    out << x;
    s = out.str();

    return s.size();
}

/**
 * Returns a truncated num/den.
 */
RT * truncateFraction(string num, string den, int size) {
    std::string num_t, den_t, x_t;
    std::stringstream out;
    RT * res;

    num_t = num.substr(0, num.size() - size);
    den_t = den.substr(0, den.size() - size);

    out << num_t << "/" << den_t;
    x_t = out.str();
    out.clear();
    out.str(std::string());

    res = new RT(x_t);

    return res;
}

/**
 * Finds a simpler point inside the polygon.
 */
Point findSimplerInteriorPoint(PolygonWithHoles polWHoles, Point p) {
    Point * p_trunc = NULL;
    int num_len, den_len, trunc_len;
    std::string num, den;
    std::stringstream out;
    RT x, y;
    RT * novo_x = NULL;
    RT * novo_y = NULL;
    Polygon pol = polWHoles.outer_boundary();

    x = p.x();
    y = p.y();

    out << x.numerator();
    num = out.str();
    out.str(std::string());
    out << x.denominator();
    den = out.str();
    out.str(std::string());

    if (!num.find_first_of("-"))
        num_len = num.size() - 1;
    else
        num_len = num.size();

    den_len = den.size();
    trunc_len = num_len < den_len ? num_len : den_len; 

    /* Searches for simplest point inside the polygon */
    for (int i=(trunc_len-1); i>=0; i--) {
        novo_x = truncateFraction(num, den, i);    
        p_trunc =  new Point(*novo_x, y);
        if (pol.bounded_side(*p_trunc) == CGAL::ON_BOUNDED_SIDE) {
            bool outOfHoles = true;

            for(PolygonWithHoles::Hole_const_iterator holeIt = polWHoles.holes_begin(); holeIt != polWHoles.holes_end(); ++holeIt) {
                if (holeIt->bounded_side(*p_trunc) != CGAL::ON_UNBOUNDED_SIDE) {
                    outOfHoles = false;
                    break;
                }
            }

            if (outOfHoles == true)
                break;
        }
            
        delete novo_x;
        delete p_trunc;
    }

    out << y.numerator();
    num = out.str();
    out.str(std::string());
    out << y.denominator();
    den = out.str();
    out.str(std::string());

    if (!num.find_first_of("-"))
        num_len = num.size() - 1;
    else
        num_len = num.size();

    den_len = den.size();
    trunc_len = num_len < den_len ? num_len : den_len; 

    for(int i=(trunc_len-1); i>=0; i--){
        delete p_trunc;
        
        novo_y = truncateFraction(num, den, i);   
        p_trunc =  new Point((*novo_x), (*novo_y));
        if (pol.bounded_side(*p_trunc) == CGAL::ON_BOUNDED_SIDE){
            bool outOfHoles = true;
            
            for(PolygonWithHoles::Hole_const_iterator holeIt = polWHoles.holes_begin(); holeIt != polWHoles.holes_end(); ++holeIt) {
                if (holeIt->bounded_side(*p_trunc) != CGAL::ON_UNBOUNDED_SIDE) {
                    outOfHoles = false;
                    break;
                }
            }
            
            if (outOfHoles == true)
                break;
        }
            
        delete novo_y;
    }

    Point res = (*p_trunc);
    delete p_trunc;
    delete novo_x;
    delete novo_y;

    return res;
}

/**
 * Finds a simpler point inside the polygon.
 */
Point findSimplerInteriorPoint(Polygon pol, Point p) {
    Point * p_trunc = NULL;
    int num_len, den_len, trunc_len;
    std::string num, den;
    std::stringstream out;
    RT x, y;
    RT * novo_x;
    RT * novo_y;

    x = p.x();
    y = p.y();

    out << x.numerator();
    num = out.str();
    out.str(std::string());
    out << x.denominator();
    den = out.str();
    out.str(std::string());

    if (!num.find_first_of("-"))
        num_len = num.size() - 1;
    else
        num_len = num.size();

    den_len = den.size();
    trunc_len = num_len < den_len ? num_len : den_len; 

    /* Searches for simplest point inside the polygon */
    for(int i=(trunc_len-1); i>=0; i--){
        novo_x = truncateFraction(num, den, i);    
        p_trunc =  new Point(*novo_x, y);
        if (pol.bounded_side(*p_trunc) == CGAL::ON_BOUNDED_SIDE)
            break;
        
        delete p_trunc;
        delete novo_x;
    }

    out << y.numerator();
    num = out.str();
    out.str(std::string());
    out << y.denominator();
    den = out.str();
    out.str(std::string());

    if (!num.find_first_of("-"))
        num_len = num.size() - 1;
    else
        num_len = num.size();

    den_len = den.size();
    trunc_len = num_len < den_len ? num_len : den_len; 

    for(int i=(trunc_len-1); i>=0; i--){
        delete p_trunc;
        
        novo_y = truncateFraction(num, den, i);    
        p_trunc =  new Point((*novo_x), (*novo_y));
        if (pol.bounded_side(*p_trunc) == CGAL::ON_BOUNDED_SIDE)
            break;
        
        delete novo_y;
    }

    Point res = (*p_trunc);
    delete p_trunc;
    delete novo_x;
    delete novo_y;
    
    return res;
}

