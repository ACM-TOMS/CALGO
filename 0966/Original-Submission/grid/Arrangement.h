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


#ifndef MYARRANGEMENT_H
#define MYARRANGEMENT_H

#include "IGrid.h"

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_observer.h>

typedef CGAL::Arr_segment_traits_2<Extended_kernel> ArrTraits;
typedef CGAL::Arr_extended_dcel<ArrTraits, bool, short, bool> DCEL; // Vertex, HE, Face
typedef CGAL::Arrangement_2<ArrTraits, DCEL> Arrangement;

/**
 * Class implemented in order to change the default behavior of the events
 * during the construction of an arrangement. With this class, it is possible
 * to save important information which is later used to verify if a cell is a
 * light or a shadow AVP. 
 */
class MyObserver : public CGAL::Arr_observer<Arrangement> {
    private:
        short _nextVal;
        short _nextValTwin;
        bool  _nextBool;
        bool  _starterEdge;
        Segment _originalSeg;

    public: 
        /**
         * Constructor. Informs that the first edges that will be added are
         * from the boundary of the polygon.
         */
        MyObserver(Arrangement& arr) : CGAL::Arr_observer<Arrangement>(arr) { _starterEdge = true; }

        /**
         * Sets variable used to inform if the next edges to be added are from
         * the boundary or induced by visibility polygons.
         */
        void setStarterEdge(bool b) { _starterEdge = b; }

        /**
         * Sets variable responsible for keeping the original edge being
         * currently included in the arrangement.
         */
        void setOriginalSeg(Segment seg) { _originalSeg = seg; }

        /**
         * Function called before creating an edge. Keeps information about 
         * the visible side of the edge.
         */
        virtual void before_create_edge (const X_monotone_curve_2& s, Vertex_handle v1, Vertex_handle v2) {
            if (_starterEdge) return;
            Segment sIni = (Segment) s;
            if (Segment(v1->point(), v2->point()).direction() != sIni.direction()) {
                _nextVal = 1;
                _nextValTwin = 0;
            } else {
                _nextVal = 0;
                _nextValTwin = 1;
            }
        }

        /**
         * Function called after creating an edge. Sets information about 
         * the visible side of the edge.
         */
        virtual void after_create_edge (Halfedge_handle e) {
            if (!_starterEdge) {
                e->set_data(_nextVal);
                e->twin()->set_data(_nextValTwin);
            } else {
                // polygon edges must have -1!!!!
                e->set_data(-1);
                e->twin()->set_data(-1);
            }
        }

        /**
         * Function called before spliting an existing edge. Saves information
         * about the visible side of the edge (which must not change).
         */
        virtual void before_split_edge (Halfedge_handle e, Vertex_handle v, const X_monotone_curve_2& s1, const X_monotone_curve_2& s2){
            _nextVal = e->data();
            _nextValTwin = e->twin()->data();
        }

        /**
         * Function called after spliting an existing edge. Sets information
         * about the visible side of both edges which were induced by the 
         * split event. The new edges should keep the same visibility status
         * of the original one.
         */
        virtual void after_split_edge (Halfedge_handle e1, Halfedge_handle e2) {
            e1->set_data(_nextVal);
            e2->set_data(_nextVal);
            e1->twin()->set_data(_nextValTwin);
            e2->twin()->set_data(_nextValTwin);
        }

        /**
         * Function called before modifying an edge. In our case, it means that
         * another edge is being placed over an existing one. It is necessary
         * to verify the direction of both edges to see if one or both sides
         * will be visible.
         */      
        virtual void before_modify_edge (Halfedge_handle e, const X_monotone_curve_2& s) {
            // especial edges (those who are light and shadow) also have -1
            _nextBool = false;
            if (Segment(e->source()->point(), e->target()->point()).direction() != _originalSeg.direction()) {
                if (e->data() == 0) {
                    _nextBool = true;
                    _nextVal = -1;
                    _nextValTwin = -1;
                }
            } else {
                if (e->data() == 1) {
                    _nextBool = true;
                    _nextVal = -1;
                    _nextValTwin = -1;
                }
            }
        }

        /**
         * Function called before modifying an edge. In our case, it means that
         * another edge is being placed over an existing one. Sets the information
         * saved before.
         */
        virtual void after_modify_edge (Halfedge_handle e) {
            if (_nextBool) {
                e->set_data(_nextVal);
                e->twin()->set_data(_nextValTwin);
            }
        }
};

#endif
