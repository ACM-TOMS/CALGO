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


#include <iostream>
#include <iomanip>
#include <string>
#include "PolygonExt.h" 
#include "polAlgorithms.h"

#define OP_FINISH    1
#define OP_ADVANCE   2
#define OP_RETARD    4
#define OP_SCAN      8

/**
 * Finds the closest point on the boundary to z.
 */
Point PolygonExt::getClosestOnBoundary(Point z) {
    Point s;
    RT distMin = -1;

    Polygon::Edge_const_iterator eIt = this->edges_begin();
    for (; eIt != this->edges_end(); ++eIt) {
        // gets the z projection on the segment eIt
        Point aux = eIt->supporting_line().projection(z);
        // verify if the projection is not on the segment, so made one of
        // the segments vertex the closest point.
        if (!eIt->has_on(aux)) {
            RT dist0 = CGAL::squared_distance(eIt->vertex(0), z);
            RT dist1 = CGAL::squared_distance(eIt->vertex(1), z);
            if (dist0 < dist1) {
                aux = eIt->vertex(0);
            } else {
                aux = eIt->vertex(1);
            }
        } 
        RT auxDist = CGAL::squared_distance(aux, z);
        if (auxDist < distMin || distMin == -1) {
            s = aux;
            distMin = auxDist;
        }
    }
    
    return s;
}


/**
 * Computes the visibility polygon of a point z in P.
 * Algorithm: B. Joe and R. B. Simpson. Visibility of a simple polygon from a point. 
 * Report CS-85-38, Dept. Math. Comput. Sci., Drexel Univ., Philadelphia, PA, 1985.
 */
PolygonExt PolygonExt::getVisibility(Point z) {

    bool debug = false;

    if (! this->is_counterclockwise_oriented()) {
        this->reverse_orientation();
    }

    std::vector<Point> vis;
    std::vector<Point> pol;
    std::vector<RT> alpha;
    std::vector<RT> alphaS;
    int oper = OP_ADVANCE;
    bool ccw = true;
    bool outNormal = true;
    int ind = 0;
    int top = 0;
    int indJ = 0;
    Point w;
    Point zOnBound = getClosestOnBoundary(z);
    PolygonExt ret; 

    RT xmin = (* this->left_vertex()).x();
    RT xmax = (* this->right_vertex()).x();
    RT ymax = (* this->top_vertex()).y();
    RT ymin = (* this->bottom_vertex()).y();
    
    Polygon pbox;
    pbox.push_back(Point(xmin - 1, ymin - 1));
    pbox.push_back(Point(xmin - 1, ymax + 1));
    pbox.push_back(Point(xmax + 1, ymax + 1));
    pbox.push_back(Point(xmax + 1, ymin - 1));
    
    if (debug) {
        std::cout << std::endl << std::endl << std::endl;
        std::cout << "POLYGON:   " << *this << std::endl;
        std::cout << "POINT:     " << z << std::endl;
    }

    if (z == zOnBound) {
        if (debug) { std::cout << ">> z == zOnBound" << std::endl; }
        // means that $z$ is on boundary of $P$. Check if it's a vertex or not and mount
        // the vertex of $P$ as $z, v_{0}, v_{1}, \cdots, v_{n}$ where $v_{n}$ is the 
        // predecessor of $z$
        Polygon::Vertex_circulator vc = this->vertices_circulator();

        bool vertex = false;
        Polygon::Vertex_iterator vIt = this->vertices_begin();
        for (; vIt != this->vertices_end(); vIt++) {
            if ((*vIt) == zOnBound) {
                vertex = true;
            break;
        }
                                                }
        if (vertex) {
            for(;; ++vc) {
                if ((*vc) == zOnBound) {
                    break;
                }
            }
            Point start = (*vc);
            for(++vc; (*vc) != start; ++vc) {
                pol.push_back(*vc);
            }
            ret.push_back(z);
        } else {
            Polygon::Edge_const_iterator ei = this->edges_begin();
            for (; ei != this->edges_end(); ++ei) {
                if (ei->has_on(zOnBound)) break;
            }
            for(;; ++vc) {
                if ((*vc) == ei->vertex(1)) break;
            }
            Point start = (*vc);
            pol.push_back(*vc);
            for(++vc; (*vc) != start; ++vc) {
                pol.push_back(*vc);
            }
        }
    } else {
        if (debug) { std::cout << ">> z != zOnBound" << std::endl; }
        // means that $z$ is in the interior of $P$. we denote the vertex of $P$ as
        // $v_{0}, v_{1}, \cdots, v_{n-1}, v_{n}$ where $v_{n} = v_{0}$
        Polygon::Vertex_circulator vc = this->vertices_circulator();
        
        bool vertex = false;
        Polygon::Vertex_iterator vIt = this->vertices_begin();
        for (; vIt != this->vertices_end(); vIt++) {
            if ((*vIt) == zOnBound) {
                vertex = true;
            break;
        }
                                                }
        if (vertex) {
            for(;; ++vc) {
                if ((*vc) == zOnBound) {
                    break;
                }
            }
        } else {
            pol.push_back(zOnBound);
            Polygon::Edge_const_iterator ei = this->edges_begin();
            for (; ei != this->edges_end(); ++ei) {
                if (ei->has_on(zOnBound)) break;
            }
            for(;; ++vc) {
                if ((*vc) == ei->vertex(1)) break;
            }
        }
         
        Point start = (*vc);
        pol.push_back(*vc);
        for(++vc; (*vc) != start; ++vc) {
            pol.push_back(*vc);
        }
        if (vertex) {
            pol.push_back(*vc);
        } else {
            pol.push_back(zOnBound);
        }
    }

    alpha.push_back(0.0);
    Point vAux = pol.at(0);
    int n = pol.size();
    n--;

    if (debug) {
        std::cout << std::endl << "ALPHA ANGLES: (" << pol.at(0) << "):" << 0.0 << " ";
        std::cout << std::endl;
    }
    
    for(int i = 1; i <= n; i++) {
        RT alphaVI = getAlpha(alpha.at(i-1), vAux, z, pol.at(i));
        alpha.push_back(alphaVI);
        if (debug) {
            std::cout << "(" << pol.at(i) << "):" << CGAL::to_double(alphaVI) << " || " << CGAL::to_double(getPseudoAngle(z,vAux)) \
                      << " - " << CGAL::to_double(getPseudoAngle(z,pol.at(i))) << " ## ";
            std::cout << std::endl;
        }
        vAux = pol.at(i);
    }
    if (debug) {
        std::cout << std::endl;
    }

    vis.push_back(pol.at(0));
    alphaS.push_back(alpha.at(0));

    if (alpha.at(1) >= alpha.at(0)) {
        oper = OP_ADVANCE;
    } else {
        oper = OP_SCAN;
        ccw = true;
        w = pointOnInfDir(z, pol.at(0), pbox);
    }

    
    for (; oper != OP_FINISH; ) {
            if (debug) {
                std::cout << std::endl << std::endl << "CUR VIS POL: ";
                for(int itVis = 0; itVis < (int) vis.size(); itVis++) {
                    std::cout << "(" << vis.at(itVis) << ") ";
                }
                std::cout << std::endl;
            }
        switch(oper) {
            case OP_ADVANCE:
                if (debug) {
                    std::cout << "OP_ADVANCE: " << ind+1 << "(" << pol.at(ind+1) << ")" << "  " << top << "(" << vis.at(top) << ")" << std::endl;
//                  std::cout << 2 * M_PI << "  " << alpha.at(ind+1) << std::endl;
                    int tmpSize = vis.size();
                    for (int i = 0; i < tmpSize; i++ ) {
                        std::cout << "(" << vis.at(i) << ") ";
                    }
                    std::cout << std::endl;
                }
                if (alpha.at(ind+1) <= RT(8)) {
                    ind++;
                    top++;
                    vis.push_back(pol.at(ind));
                    alphaS.push_back(alpha.at(ind)); 
                    if (ind == n) {
                        oper = OP_FINISH;
                    } else if (alpha.at(ind+1) < alpha.at(ind) && CGAL::right_turn(pol.at(ind-1), pol.at(ind), pol.at(ind+1))) {
                        oper = OP_SCAN;
                        ccw = true;
                        w = pointOnInfDir(z, pol.at(ind), pbox);
                    } else if (alpha.at(ind+1) < alpha.at(ind) && CGAL::left_turn(pol.at(ind-1), pol.at(ind), pol.at(ind+1))) {
                        oper = OP_RETARD;
                    }
                } else {
                    if (alphaS.at(top) < RT(8)) {
                        CGAL::Object obj = CGAL::intersection(Segment(pol.at(ind), pol.at(ind+1)), Ray(z, pol.at(0)));
                        if (const Point *point = CGAL::object_cast<Point>(&obj)) {
                            top++;
                            vis.push_back(*point);
                            alphaS.push_back(RT(8)); 
                        }
                    }
                    oper = OP_SCAN;
                    ccw = false;
                    w = pol.at(0);
                }
                break;
            case OP_RETARD:
                if (debug) {
                    std::cout << "OP_RETARD: " << ind+1 << "(" << pol.at(ind+1) << ")" << "  " << top << "(" << vis.at(top) << ")" << std::endl;
                }
                outNormal = false;
                for (indJ = top - 1; indJ >= 0; indJ--) { // antes era indJ > 0
                    if (debug) {
                        std::cout << indJ << " " << vis.at(indJ) << ":" << alphaS.at(indJ) << " " << pol.at(ind+1) << ":" << alpha.at(ind+1) << " " << vis.at(indJ+1) << ":" << alphaS.at(indJ+1) << std::endl;
                        std::cout << indJ << " " << alpha.at(ind+1) << " " << alphaS.at(indJ) << "  " << alphaS.at(indJ+1) << std::endl;
                    }
                    if (alphaS.at(indJ) < alpha.at(ind+1) && alpha.at(ind+1) <= alphaS.at(indJ+1)) {
                        if (debug) {
                            std::cout << "CASE RETARD - A" << std::endl;
                        }
                        outNormal = true;
                        break;
                    }
                    if (alpha.at(ind+1) <= alphaS.at(indJ) && alphaS.at(indJ) == alphaS.at(indJ+1) && pol.at(ind) != vis.at(indJ+1)) {
                        CGAL::Object obj = CGAL::intersection(Segment(pol.at(ind), pol.at(ind+1)), Segment(vis.at(indJ), vis.at(indJ+1)));
                        if (const Point *point = CGAL::object_cast<Point>(&obj)) {
                            w = (*point);
                            if (debug) {
                                std::cout << "CASE RETARD - B, point w equal to " << w << std::endl;
                            }
                            outNormal = true;
                            break;
                        }
                    }
                    vis.pop_back();
                    alphaS.pop_back();
                }
                top = indJ;
                if (debug) {
                    std::cout << "INDJ:  " << indJ << std::endl;
                }
                if (alphaS.at(top) < alpha.at(ind+1)) {
                    ind++;
                    CGAL::Object obj = CGAL::intersection(Segment(vis.at(top), vis.at(top+1)), Ray(z, pol.at(ind)));
                    vis.pop_back();
                    alphaS.pop_back();
                    if (const Point *point = CGAL::object_cast<Point>(&obj)) {
                        top++;
                        vis.push_back(*point);
                        alphaS.push_back(alpha.at(ind));
//                      alphaS.push_back(getAlpha(alphaS.at(top-1), vis.at(top-1), z, *point));
                    }
                    top++;
                    vis.push_back(pol.at(ind));
                    alphaS.push_back(alpha.at(ind));
                    if (ind == n) {
                        oper = OP_FINISH;
                    } else if (alpha.at(ind+1) >= alpha.at(ind) &&  CGAL::right_turn(pol.at(ind-1), pol.at(ind), pol.at(ind+1))) {
                        oper = OP_ADVANCE;
                    } else if (alpha.at(ind+1) > alpha.at(ind) && CGAL::left_turn(pol.at(ind-1), pol.at(ind), pol.at(ind+1))) {
                        oper = OP_SCAN;
                        ccw = false;
                        w = pol.at(ind);
                        top--;
                        vis.pop_back();
                        alphaS.pop_back();
                    } else {
                        top--;
                        vis.pop_back();
                        alphaS.pop_back();
                    }
                } else {
                    vis.pop_back();
                    alphaS.pop_back();
                    if (debug) {
                        std::cout << "RETARD - ELSE CASE B : " << pol.at(ind) << ":" << alpha.at(ind) << " " << pol.at(ind+1) << ":" << alpha.at(ind+1) << " " << pol.at(ind+2)<< ":" << alpha.at(ind+2) << std::endl;
                    }
                    if (alpha.at(ind+1) == alphaS.at(top) && alpha.at(ind+2) > alpha.at(ind+1) && CGAL::right_turn(pol.at(ind), pol.at(ind+1), pol.at(ind+2))) {
                        oper = OP_ADVANCE;
                        ind++;
                        top++;
                        vis.push_back(pol.at(ind));
                        alphaS.push_back(alpha.at(ind));
                    } else {
                        oper = OP_SCAN;
                        ccw = true;
                        if (!outNormal) {
                            outNormal = true;
                            w = pointOnInfDir(z, vis.at(top), pbox);
                        }
                    }
                }
                break;
            case OP_SCAN:
                ind++;
                if (ind == n) {
                    oper = OP_FINISH; // to avoid out_of_range
                    break;
                }
                if (debug) {
                    std::cout << "OP_SCAN: " << ind+1 << "(" << pol.at(ind+1) << ")" << "  " << top << "(" << vis.at(top) << ")" << std::endl;
                }
                if (ccw && alpha.at(ind+1) > alphaS.at(top) && alphaS.at(top) >= alpha.at(ind)) {
                    if (debug) {
                        std::cout << "SCAN CCW " << ind << " "  << vis.size() << " " << top << " " << pol.at(ind+1) << ":" << alpha.at(ind+1) << " " << vis.at(top) << ":" << alphaS.at(top) << " " << pol.at(ind) << ":" << alpha.at(ind) << " " << w << std::endl;
                    }
                    CGAL::Object obj = CGAL::intersection(Segment(pol.at(ind), pol.at(ind+1)), Segment(vis.at(top), w));
                    if (const Point *point = CGAL::object_cast<Point>(&obj)) {
                        if ((alpha.at(ind-1) <= alpha.at(ind) && *point != pol.at(ind+1)) || 
                            (alpha.at(ind-1) > alpha.at(ind) && *point != pol.at(ind))) {
                            if (debug) {
                                std::cout << "SCAN CCW IF " << std::endl;
                            }
                            top++;
                            vis.push_back(*point);
                            alphaS.push_back(alphaS.at(top-1));
                            oper = OP_ADVANCE;
                        }
                    }
                } else {
                    if (debug) {
                        std::cout << ind << " "  << vis.size() << " " << top << " " << pol.at(ind+1) << ":" << alpha.at(ind+1) << " " << vis.at(top) << ":" << alphaS.at(top) << " " << pol.at(ind) << ":" << alpha.at(ind) << " CCW: " << ccw << std::endl;
                    }
                    if(!ccw && alpha.at(ind+1) <= alphaS.at(top) && alphaS.at(top) < alpha.at(ind)) {
                        CGAL::Object obj = CGAL::intersection(Segment(pol.at(ind), pol.at(ind+1)), Segment(vis.at(top), w));
                        if (CGAL::object_cast<Point>(&obj)) {
                            oper = OP_RETARD;
                        }
                    }
                }
                break;
            default:
                break;
        }
    }

    int visN = vis.size();

    for (int i = 0; i < visN; i++) {
        ret.push_back(vis.at(i));
    }

    return cleanPol(ret);
}

